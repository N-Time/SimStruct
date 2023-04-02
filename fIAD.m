% Interlayer (Interstorey) Area Damage, 层间面积损伤
% UPDATE: Compared to 'dm_fIAD.m',
% 1) implemented by the MATLAB build-in function 'zerocrossrate' (Introduced in R2021b);
% 2) consider the negative response of interlayer relative deformation,
% indicated by the negative value LD.

% %%%% LOADING CASE INPUT
% loadMainPath = 'E:\ANSYS\MAS_IDA\LA_FV_no_Gravity_mat\LF10NGpr'; %
% % database file main path
% databaseName = 'LF10NGpr';  % database name
% component0 = 'CALayer';   % component name
% gm_dirc = 'Y';   % GM incidence
% dis_dirc = 'Y';   % Response (displacement) direction
% disp_cat = 'Peak';   % database label: Peak or Residual 
% 
% num_record = 37;  % record num. in the database
% imLevel = 1.8;   % IM level
% dIM = 0.2;   % IM level step
% pl_sw = 1;   % plot switch: 1 = on, others = off
% LSdataFile = "D:\Wen\Matlab\Project\ModularizedDamageModel\IADM\LS_RBS.mat";  % LS table data
% 
% % Retrieve xND at specific (L,d,SP,zMB,t) in DLd(SP,xND,zMB,t)
% EDPdata = load_DMdata_LF10NGpr(loadMainPath, gm_dirc, databaseName,...
%     component0, disp_cat);
% [xx, yy] = LF10NGpr_disp_span(EDPdata,component0,dis_dirc,...
%                 num_record,imLevel,dIM,0,LSdataFile);
% 
% %%%% IAD
% dirc_id = ANSYS_dirc_label(dis_dirc);  % direction id in ANSYS
% [num_node,~,~] = ANSYS_com_num_node(component0,dirc_id,5,51); % node total num. for the given component
% miu_s = 0.25;  %%%%% KEY parameter of the onset of the interlayer failure
% LS_RBS = load(LSdataFile);  % predetermined LSs table
% Limit_State = LS_RBS.Limit_State;
% len_ele = 0.645;   % element length in FEM
% 
% % Calculate IAD result with the input of  xND at specific (L,d,SP,zMB,t) in DLd(SP,xND,zMB,t)
% yy1 = ([yy(1:125),-yy(126:end)]); % NxM delta(x), N is the samples, M is the input channels (response components)
% x = [yy', yy1'];
% d = xx'; % the location coordinate of x
% signalLevel = Limit_State{1,component0};  % threshold by yield disp.
% rspVarName = ["ResponseX","ResponseY"];
% 
% % T = fIAD(d,x,signalLevel,"rspVarName",rspVarName,"isPlot",1);
% T = fIAD(d,x,signalLevel,"rspVarName",rspVarName);

function T = fIAD(d,x,signalLevel,varargin)

    % Default parameters
    ip = inputParser;   % Parser
    % default
    addParameter(ip,'rspVarName',("comp" + (1:size(x,2)))');  % the row name of DI table
    addParameter(ip,'isPlot',0);  % the time var name in 'dr'
    parse(ip,varargin{:});  % update
    % read the par.
    rspVarName = ip.Results.rspVarName;
    isPlot = ip.Results.isPlot;
    
    % Positive IAD
    [T_ps, sc_ps] = fIAD0(d,x,signalLevel,"rspVarName",rspVarName);
    % Negative IAD
    [T_ng, sc_ng] = fIAD0(d,x,-signalLevel,"rspVarName",rspVarName);

    % Merged + - IAD
    rspVarName = rspVarName';
    xL = cell(size(x,2),1);
    xR = cell(size(x,2),1);
    xloc = cell(size(x,2),1);
    LD = cell(size(x,2),1);
    for i = 1:1:size(x,2)
        xL{i} = [T_ps.xL{i}; T_ng.xL{i}];
        xR{i} = [T_ps.xR{i}; T_ng.xR{i}];
        xloc{i} = [T_ps.xloc{i}; T_ng.xloc{i}];
        LD{i} = [T_ps.LD{i}; T_ng.LD{i}];
    end
    T = table(rspVarName, xL, xR, xloc, LD);


    % Plot
    if isPlot == 1
        figure

        plot(d,x)
        hold on

        scatter([sc_ps.BC; sc_ng.BC], [sc_ps.XX; sc_ng.XX],'o','filled','r')
        
        ax = gca;
        xlabel('Span / m');
        ylabel('\it {\delta(x)}\rm m');
        set(gca,'fontsize',18,'fontname','Times');
%         set(gca,'linewidth',2)

        xlim([0, d(end)]) % limit the span
%         ylim([-3, 3])

        grid on
        box on
    end


end


function [T0,sc] = fIAD0(d,x,signalLevel,varargin)
    % by 'zerocrossrate': (b) ii), iii), and (c) i), ii), iii) for conservative boundaries
    % 'zerocrossrate' returns N*M matrix, N = num of windows, M = num of input channels
    % (b) ii), iii) consider 0 as positive
    % (c) i) left BC = pick up the left point; ii) right BC = pick up the right point
    % (c) iii) add the right boundary

    
    d = repmat(d,1,size(x,2));
    
    % The zerocrossrate function returns an index at the 'next' sample following a crossing.
    % (c) rude rising
    [~,~,idRise] = zerocrossrate(x, ...
        "Level",signalLevel,"ZeroPositive",1,"TransitionEdge","rising");  % only rising
    if size(size(idRise),2) > 2
        idRiseM0 = squeeze(idRise);  % multi-channel input => 1xNxM
    else
        idRiseM0 = idRise';  % single-channel input => 1xN
    end
    % (c) rude falling
    [~,~,idFall] = zerocrossrate(x, ...
        "Level",signalLevel,"ZeroPositive",1,"TransitionEdge","falling");  % only falling
    if size(size(idFall),2) > 2
        idFallM0 = squeeze(idFall);  % multi-channel input => 1xNxM
    else
        idFallM0 = idFall';  % single-channel input => 1xN
    end

    % (c) Conservative Criterion: i) left BC and ii) right BC
    if signalLevel >= 0
        % rising = left BC: id - 1; falling = right BC: auto to pin at the 'next' sample index
        idRiseM1 = [idRiseM0(2:end,:); false(1,size(idRiseM0,2))];  % move forward 1 element and add 0 at the last row
        idRiseM1(1,:) = idRiseM0(1,:) | idRiseM1(1,:); % union the 1st row and the original 1st row

        idFallM1 = idFallM0;
    else  % signalLevel < 0
        % rising = right BC: auto to pin at the 'next' sample index; falling = left BC: id - 1
        idRiseM1 = idRiseM0;  

        idFallM1 = [idFallM0(2:end,:); false(1,size(idFallM0,2))];  % move back 1 element and add 0 at the 1st row
        idFallM1(1,:) = idFallM0(1,:) | idFallM1(1,:); % union the last row as the original last row
    end
    
    
    % (c) iii) add the right boundary
    if signalLevel >= 0  % rising = left BC; falling = right BC
        for ii = 1:1:size(x,2)  % for each column (input channel)
            xCol = x(:,ii);
            if size(xCol(idRiseM1(:,ii)),1) > size(xCol(idFallM1(:,ii)),1)  % for the case without the right BC
                idFallM1(end,ii) = true;  % add the index of the supplement right BC
            end
        end

    else  % signalLevel < 0  % rising = right BC; falling = left BC
        for ii = 1:1:size(x,2)  % for each column (input channel)
            xCol = x(:,ii);
            if size(xCol(idFallM1(:,ii)),1) > size(xCol(idRiseM1(:,ii)),1)  % for the case without the right BC
                idRiseM1(end,ii) = true;  % add the index of the supplement right BC
            end
        end
    end


    
%     % idBoth = idRiseM1 | idFallM1;
%     idBoth = idRiseM1 | idFallM0;

    % Output
    results = struct;
    xL = cell(1,size(x,2));  % for each channel
    xR = cell(1,size(x,2));
    xloc = cell(1,size(x,2));
    LD = cell(1,size(x,2));
    for ii = 1:1:size(x,2)  % for each channel
        %%% (*d*) UPDATE: select the true 'rise' or 'fall' to the sign of response
        if signalLevel >= 0
            xL{ii} = d(idRiseM1(:,ii),ii);
            xR{ii} = d(idFallM1(:,ii),ii);
        else
            xL{ii} = d(idFallM1(:,ii),ii);
            xR{ii} = d(idRiseM1(:,ii),ii);
        end

        % Left BC
        if isempty(xL{ii})
            xL{ii} = [];  % better to table operation
        end

        % Right BC
        if isempty(xR{ii})
            xR{ii} = [];
        end
        
        if isempty(xL{ii}) || isempty(xR{ii})
            xloc{ii} = [];
            LD{ii} = [];
        else
            
            %%% (*d*) UPDATE: determine the sign of the response (+ or -)
            if size(xL{ii},1) == size(xR{ii},1)
                LD{ii} = d(idFallM1(:,ii),ii) - d(idRiseM1(:,ii),ii);
                xloc{ii} = (xR{ii} - xL{ii})./2 + xL{ii};
            else
                error("An inconsistent left and right BCs are determined at Response Channel " + ii)
            end

        end
    end

%     results.xL = xL;
%     results.xR = xR;
%     results.xloc = xloc;
%     results.LD = LD;
    
    
    T0.xL = xL';
    T0.xR = xR';
    T0.xloc = xloc';
    T0.LD = LD';
%     T = table(xL, xR, LD, xloc, 'RowNames',rspVarName);
%     rspVarName = rspVarName';
%     T0 = table(rspVarName, xL, xR, xloc, LD);


    sc.BC = d(idRiseM1 | idFallM1);
    sc.XX = x(idRiseM1 | idFallM1);
%     % Plot
%     if isPlot == 1
%         figure
% 
%         plot(d,x)
%         hold on
% 
%         scatter(d(idRiseM1), x(idRiseM1),'<','filled','r')
%         hold on
%         scatter(d(idFallM1), x(idFallM1),'>','filled','r')
%         
%         ax = gca;
%         xlabel('Span / m');
%         ylabel('\it {\delta(x)}\rm m');
%         set(gca,'fontsize',18,'fontname','Times');
% %         set(gca,'linewidth',2)
% 
%         xlim([0, d(end)]) % limit the span
% %         ylim([-3, 3])
% 
%         grid on
%         box on
%     end

end