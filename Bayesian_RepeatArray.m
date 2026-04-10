function C=Bayesian_RepeatArray(B, num)
if ~isvector(B)
    error('输入必须是一维数组');
end
C = repmat(B, 1, num);
end