function R = white_metric(X,m)
%WHITE_METRIC	 Compute a metric for comparing the relative "whiteness" of
%                a random process.
%   R = WHITE_METRIC(X,m) computes and returns a whiteness metric for the
%   random sequence X based on n elements of the autocovariance sequence
%   for X
%
%   See The Digital Signal Processing Handbook 16.3.2

    [numx, lenx] = size(X);
%     if (lenx > numx)
%         X = X';
%         n = numx;
%         numx = lenx;
%         lenx = n;
%     end
    
    for r = numx:-1:1
        r_hat(r,:) = xcorr(X(r,:));
        lenr = length(r_hat(r,:));
        lag = ceil(lenr/2)+1 + [1:m];
        R(r,:) = lenx * sum(r_hat(r,lag).^2) / r_hat(r,ceil(lenr/2)+1).^2;
    end
    
    
end