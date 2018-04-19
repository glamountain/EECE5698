function [cxcor, lags] = xcorr_circ(a,b)

% Normalize
a = a(:) / norm(a);
b = b(:) / norm(b);

cxcor = zeros(1,length(b));
% Shift and compute cross correlation at each shift
for k=1:length(b)
    
    cxcor(k) = a' * b;
    b = circshift(b,1,1);

end
lags = [0:length(b)-1];