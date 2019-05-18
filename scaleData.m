function [Bs,ms] = scaleData(B,m)
% normalize data to [-1,1]

ub = max([B(:); m(:)]);
lb = min([B(:); m(:)]);

% scale to [-1, 1]
Bs = ((B-lb)./ub)*2-1;
ms = ((m-lb)./ub)*2-1;

% scale to [0, 1]
% Bs = ((B-lb)./ub);
% ms = ((m-lb)./ub);

