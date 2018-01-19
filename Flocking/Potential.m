% A helper script that takes a combination of radii as input and returns a
% potential function value
% Author: Logan Beaver
% Updated: 1/19/2018
function [V, dv] = Potential(r, Rhat, R)
%desired radius is Rhat
%maximum range is R
%1/x from 0 to Rhat

% minVal = 1/Rhat;
% constant = 5;
% 
% 
% if r < Rhat
%     V = 0.01/r;
% elseif r < R
%     V = (minVal - constant)/(Rhat - R)*(r-Rhat) + minVal;
% else
%     V = constant;
% end

if r <= 5
    V = 1./(r.^2) + log(r.^2);
    dv = 2*(r^2-1)/r^3;
else
    V = 1./(5^2) + log(5^2) * ones(size(r));
    dv = 0;
end

end

