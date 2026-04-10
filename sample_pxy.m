% Joint PDF
%                Also using inverse transform sampling
function [sampleX,sampleY] = sample_pxy (pxy,x,y)
% sample latitude
% Sum the latitudes for a given longitude and sample it
Fx = cumsum(sum(pxy,1));
sampleX = interpola(Fx,x,rand);
indxXUpp = find(x>sampleX,1,'first');
indxXLow = find(x<sampleX,1,'last');

if isempty(indxXLow)
    py_x = pxy(:,1);  % left bound
elseif isempty(indxXUpp)
    py_x = pxy(:,end);% right bound
else
    % Continuous spatial sampling is achieved through interp1, 
    % breaking through the limitation of grid discretization.
    py_x = interp1([x(indxXLow);x(indxXUpp)],[pxy(:,indxXLow)';pxy(:,indxXUpp)'], ...
        sampleX);
end
py_x = py_x/sum(py_x);

% sample latitude
Fy_x = cumsum(py_x);    
sampleY = interpola(Fy_x,y,rand);

end
        