% SimSturct: EDP modular: SEDM
% 以位移表征的离散分布构件的EDP变换函数

%%% Input:
% x = 构件特征响应，如桥墩-墩顶漂移，支座-水平变形
% x = [time, response]，行-时间；列-响应分量

% normMetric = 无量纲化系数，如桥墩-墩高，支座-1（1表示无需无量纲化）
% normMetric = [1]，标量

%%% Output:
% y = 最终输出的EDP
% y = [1, response], 行-1（时序最大绝对值）；列-响应分量

%%% e.g.
% x = [1, -5; 3, 3; 2, 4];  % 时间点3个，2个响应分量
% normMetric = 2;
% y = fedpSEDM(x, normMetric);  % 默认输出绝对值 outValue = "abs"
% % y = [1.5, 2.5];
% y = fedpSEDM(x, normMetric, outValue = "origin");  % 可选择输出原值
% % y = [1.5, -2.5];


function y = fedpSEDM(x, normMetric, options)

arguments
    x {mustBeNumeric}
    normMetric (1,1) {mustBeNumeric}
    options.outValue {mustBeText} = "abs"
end

xNorm = x./normMetric;  % 无量纲化

if options.outValue == "abs"
    % 输出绝对值
    y = max(abs(xNorm),[],1);  % 响应各分量的时序最大值绝对值

elseif options.outValue == "origin"
    % 输出原值
    [~, yID] = max(abs(xNorm),[],1,'linear');
    y = xNorm(yID);

else
    error("The outValue must be 'abs' or 'origin'.")
end



end