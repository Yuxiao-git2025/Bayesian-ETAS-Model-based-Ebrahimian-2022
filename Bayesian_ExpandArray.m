function A2=Bayesian_ExpandArray(A, num)
if ~isvector(A)
    error('输入必须是一维数组');
end
A = A(:)';
A2 = repelem(A, num);
end