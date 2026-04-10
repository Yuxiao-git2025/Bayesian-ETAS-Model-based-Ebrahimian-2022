function [cp1, cp2] = calculate_CI2(data, C1, C2)
if nargin < 3 || isempty(C2)
    C2 = 84; % 默认上百分位
end

if nargin < 2 || isempty(C1)
    C1 = 16; % 默认下百分位
end
% 计算百分位数
cp1 = prctile(data, C1);
cp2 = prctile(data, C2);

% 对于矩阵输入，按列计算
if ~isvector(data)
    cp1 = prctile(data, C1, 1);
    cp2 = prctile(data, C2, 1);
end

end