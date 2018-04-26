function R_est = Myers_cov(y,HPH)

L = size(y,2);
ny= size(y,1);

y_mean = mean(y,2);

if ny==1,
    R_est = 1/(L-1)*sum( (y-y_mean).^2 - (L-1)/L*HPH );
else
    R_est = zeros(ny,ny);
    for kk=1:L
        R_est = R_est + (y(:,kk)-y_mean)*(y(:,kk)-y_mean).' - (L-1)/L*HPH(:,:,kk);
    end
    R_est = R_est/(L-1);
end
