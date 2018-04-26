function [mu_est, R_est, mu_posterior, kappa_posterior, nu_posterior, Psi_posterior] = BayesianCovEst(y, mu_prior, kappa_prior, nu_prior, Psi_prior)

K  = size(y,2);
ny = size(y,1);

% aux computations
y_mean = mean(y,2);
% Psi_N = sum((y-y_mean)*(y-y_mean).');
Psi_N = zeros(ny,ny);
for kk=1:K
    Psi_N = Psi_N + (y(:,kk)-y_mean)*(y(:,kk)-y_mean).';
end

% update posterior parameters
mu_posterior = (kappa_prior*mu_prior + K*y_mean)/(kappa_prior + K);
kappa_posterior = kappa_prior + K;
nu_posterior = nu_prior + K;
Psi_posterior = Psi_prior + Psi_N + (kappa_prior*K)/(kappa_prior + K)*(y_mean - mu_prior)*(y_mean - mu_prior).' ;

% compute estimates
mu_est = mu_posterior;
if nu_posterior - ny - 1 > 0
    R_est = Psi_posterior/(nu_posterior - ny - 1);   % mean
else
    R_est = Psi_posterior/(nu_posterior + ny + 1);   % mode
end

end