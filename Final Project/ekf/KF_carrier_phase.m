%function [x_kf_est, x] = KF_carrier_phase(Niter, N)
clear
close all

Niter = 10;

% signal parameters
if ~exist('Niter', 'var')
    Niter = 10;     % Monte Carlo trials
end
if ~exist('N', 'var')
    N=3e2;         % number of samples
end

% nx = 4;
nx = 3;
ny = 2;
x = zeros(nx,N);
y = zeros(ny,N);

online = true;
Ntransient = 1e2;

T = 1;

% State Equation
f  = @(x)([1,2*pi*T,pi*T^2/2; 0,1,T; 0,0,1]*x);
Ff = @(x)([1,2*pi*T,pi*T^2/2; 0,1,T; 0,0,1]);
B  = [1,0,0; 0,1,0; 0,0,1];

% f  = @(x)([1,0,T,0; 0,1,0,T; 0,0,1,0; 0,0,0,1]*x);
% Ff = @(x)([1,0,T,0; 0,1,0,T; 0,0,1,0; 0,0,0,1]);
% B  = [T^2/2,0; 0,T^2/2; T,0; 0,T];

% Measurement Equation
a = 10;
h  = @(x)(a*[ cos(x(1,:)); sin(x(1,:))]);
Hh = @(x)(a*[-sin(x(1,:)), 0, 0; cos(x(1,:)), 0, 0]);

% h  = @(x)([1,0,0,0; 0,1,0,0]*x);
% Hh = @(x)([1,0,0,0; 0,1,0,0]);

% Error Statistics
sigma_v = 0.001;            % state standard deviation
Q       = B*sigma_v*B.';
sigma_n = 0.1;
R       = sigma_n * eye(2); % measurement covariance (actual)
% R       = diag([10,20]);


% KF parameters
x0_init_est = zeros(nx,1);
Cov0_init   = diag([1,1,1]);
%                 x0_init_est = [100;100;0;0];
%                 Cov0_init = diag([20,20,0.1,0.1]);
Q_kf_0      = Q;
% R_kf_0      = R; % 50*eye(ny);      % measurement covariance (kalman)
R_kf_0      = 10*eye(ny);
% R_kf_0 = 50*eye(ny);

% Initial Values
x(:,1) = zeros(nx,1);
y(:,1) = h(x(:,1)) + sqrt(R)*randn(ny,1);

for ii=2:N
    x(:,ii) = f(x(:,ii-1)) + B*sqrt(sigma_v)*randn(size(B,2),1);
    y(:,ii) = h(x(:,ii)) + sqrt(R)*randn(ny,1);
end

% return

% method parameters (Myers)
L = 100;        % window length

% method parameters (Bayesian Covariance Estimation - BCE)
mu_prior_0 = zeros(ny,1);
kappa_prior_0 = 1;
nu_prior_0 = 10;
Psi_prior_0 = R*(nu_prior_0 + ny + 1);  % mode around R
% Psi_prior_0 = R_kf_0*(nu_prior_0 + ny + 1);  % mode around R_kf_0

% mu_prior_0 = zeros(ny,1);
% kappa_prior_0 = 0;        % noninf: 0
% nu_prior_0 = 0;               % noninf: 0 (or -1)
% Psi_prior_0 = 1e2*eye(ny); %R_kf_0*(nu_prior_0_ref + ny + 1); %0;%R_kf_0*(nu_prior_0 + ny + 1);  % mode around R


% % reference analysis, non-informative
% mu_prior_0_ref = zeros(ny,1);
% kappa_prior_0_ref = 0;        % noninf: 0
% nu_prior_0_ref = 0;               % noninf: 0 (or -1)
% Psi_prior_0_ref = 1e2*eye(ny); %R_kf_0*(nu_prior_0_ref + ny + 1); %0;%R_kf_0*(nu_prior_0 + ny + 1);  % mode around R


% init
mu_est = zeros(ny,N);
% S_est = zeros(ny,ny,N);
R_est = zeros(ny,ny,N);
x_kf_est = zeros(size(x));
P_kf_est = zeros(nx,nx,N);
HPH = zeros(ny,ny,N);
xinnovations = zeros(size(x));
yinnovations = zeros(size(y));
R_MSE  = zeros(ny,N);
whiteness = 0;

R_est3 = zeros(ny,ny,N);
x_kf_est3 = zeros(size(x));
P_kf_est3 = zeros(nx,nx,N);
HPH3 = zeros(ny,ny,N);
xinnovations3 = zeros(size(x));
yinnovations3 = zeros(size(y));
R_MSE3  = zeros(ny,N);

for mm = 1:Niter
    mm
    
    % init
    x_est_ant = x0_init_est;
    P_est_ant = Cov0_init;
    R_kf = R_kf_0;
    Q_kf = Q_kf_0;
    
    x_est_ant3 = x0_init_est;
    P_est_ant3 = Cov0_init;
    R_kf3 = R_kf_0;
    Q_kf3 = Q_kf_0;
    
    mu_prior = mu_prior_0;
    kappa_prior = kappa_prior_0;        % noninf: 0
    nu_prior = nu_prior_0;               % noninf: -1 or 0
    Psi_prior = Psi_prior_0;  % mode around R
    %     mu_prior_ref = mu_prior_0_ref;
    %     kappa_prior_ref = kappa_prior_0_ref;        % noninf: 0
    %     nu_prior_ref = nu_prior_0_ref;               % noninf: -1 or 0
    %     Psi_prior_ref = Psi_prior_0_ref;  % mode around R
    
    for ii = 1:N
        
        % Prediction Step
        % Fk = Ff(x_est_ant);
        [x_kf_pred,  P_kf_pred]             = EKF_pred(f,Ff,Q_kf,x_est_ant,P_est_ant);
        [yinnovations(:,ii),  HPH(:,:,ii)]  = EKF_innovate(y(:,ii),h,Hh,x_kf_pred,P_kf_pred);
        
        % Fk3 = Ff(x_est_ant3);
        [x_kf_pred3, P_kf_pred3]            = EKF_pred(f,Ff,Q_kf,x_est_ant3,P_est_ant3);
        [yinnovations3(:,ii), HPH3(:,:,ii)] = EKF_innovate(y(:,ii),h,Hh,x_kf_pred3,P_kf_pred3);
        
        if (ii > Ntransient)
            %%  Bayesian Covariance Estimation
            [mu_est(:,ii), R_est(:,:,ii), mu_posterior, kappa_posterior, nu_posterior, Psi_posterior] = ...
                BayesianCovEst(yinnovations(:,ii), mu_prior, kappa_prior, nu_prior, Psi_prior);
            
            mu_prior = mu_posterior;
            kappa_prior = kappa_posterior;
            nu_prior = nu_posterior;
            Psi_prior = Psi_posterior;
            
            % TODO: Why are these next two operations in this order
            
            % update R estimate for KF
            if (online) && any(any(isinf(R_est(:,:,ii))==0))
                R_kf = R_est(:,:,ii);
            end
            
            R_est(:,:,ii) = R_est(:,:,ii) - HPH(:,:,ii);
            % R_est(:,:,ii) = R_est(:,:,ii) + yinnovations(:,ii)*yinnovations(:,ii).' - HPH(:,:,ii);
            
            %%  Myers Covariance Estimation
            if (ii>L+Ntransient)
                R_est3(:,:,ii) = Myers_cov(yinnovations3(:,ii-L+1:ii),HPH3(:,:,ii-L+1:ii));
                
                if (online)
                    R_kf3 = R_est3(:,:,ii);
                end
            else
                R_est3(:,:,ii) = R_kf3;
            end
        else
            
            R_est(:,:,ii) = R_kf;
            R_est3(:,:,ii) = R_kf3;
        end
        
        % Update Step
        Hk  = Hh(x_kf_pred);
        Hk3 = Hh(x_kf_pred3);
        [x_kf_est(:,ii), P_kf_est(:,:,ii), xinnovations(:,ii)] = KF_upd(yinnovations(:,ii),Hk,R_kf,x_kf_pred,P_kf_pred);
        [x_kf_est3(:,ii), P_kf_est3(:,:,ii), xinnovations3(:,ii)] = KF_upd(yinnovations3(:,ii),Hk3,R_kf3,x_kf_pred3,P_kf_pred3);
        
        % Note: Can use KF_upd because nonlinear residuals are computed in EKF_pred
        
        x_est_ant = x_kf_est(:,ii);
        P_est_ant = P_kf_est(:,:,ii);
        
        x_est_ant3 = x_kf_est3(:,ii);
        P_est_ant3 = P_kf_est3(:,:,ii);
    end
    
    %     if (~cov_est)
    %         % If Offline, perform Batch Estimation
    %         [mu_est_bat, R_est_bat, mu_posterior_bat, kappa_posterior_bat, nu_posterior_bat, Psi_posterior_bat] = ...
    %             BayesianCovEst(yinnovations, mu_prior_0, kappa_prior_0, nu_prior_0, Psi_prior_0);
    %     end
    
    % Check Whiteness of Innovations
    whiteness = whiteness + white_metric(yinnovations,1);
    
    % for RMSE calculation
    for dd = 1:ny
        R_MSE(dd,:)  = R_MSE(dd,:)  + (squeeze(R_est(dd,dd,:) - R(dd,dd)).^2).';
        %     R_MSE2(1,:)  = R_MSE2(1,:)  + (squeeze(R_est2(1,1,:) - R(1,1)).^2).';
        R_MSE3(dd,:)  = R_MSE3(dd,:)  + (squeeze(R_est3(dd,dd,:) - R(dd,dd)).^2).';
    end
end

whiteness = whiteness/Niter;
R_RMSE = sqrt(R_MSE/Niter);
% R_RMSE2 = sqrt(R_MSE2/Niter);
R_RMSE3 = sqrt(R_MSE3/Niter);

% figure,
% %subplot(2,1,1)
% plot3(x(1,:),x(2,:),x(3,:),'LineWidth',3), hold on,
% plot3(x_kf_est(1,:),x_kf_est(2,:),x_kf_est(3,:),'r','LineWidth',2)
% plot3(x_kf_est3(1,:),x_kf_est3(2,:),x_kf_est3(3,:),'g')
% xlabel('$\theta_{d,k}$ (rad), Reciever Dynamic Phase Variation','Interpreter','latex');
% ylabel('$f_{d,k}$ (Hz), Carrier Doppler frequency shift','Interpreter','latex');
% zlabel('$f_{r,k}$ (Hz/s), Carrier Doppler frequency rate','Interpreter','latex'), grid on
% title('State Evolution, Kalman Estimate, perfectly known state covariance');
% legend('True State Evolution','EKF State Estimate (Bayes)','EKF State Estimate (Myers)')


figure,
labels  = {'$\theta_{d,k}$ (rad), Reciever Dynamic Phase Variation', ...
    '$f_{d,k}$ (Hz), Carrier Doppler frequency shift' ...
    '$f_{r,k}$ (Hz/s), Carrier Doppler frequency rate'};
for n = 1:3
    subplot(3,1,n)
    plot(1:N,x(n,:),'LineWidth',3), hold on,
    plot(1:N,x_kf_est(n,:),'r','LineWidth',2)
    plot(1:N,x_kf_est3(n,:),'g')
    %xlabel('$\theta_{d,k}$ (rad), Reciever Dynamic Phase Variation','Interpreter','latex');
    %ylabel('$f_{d,k}$ (Hz), Carrier Doppler frequency shift','Interpreter','latex'), grid on
    title(labels{n},'Interpreter','latex');
    legend('True State Evolution','EKF State Estimate (Bayes)','EKF State Estimate (Myers)')
end

% %subplot(2,1,2)
% figure,
% y_est = h(x_kf_est);
% y_est3 = h(x_kf_est3);
% plot(y(1,:),y(2,:),'*g'), hold on, plot(y_est(1,:),y_est(2,:),'ok'),hold on, plot(y_est3(1,:),y_est3(2,:),'ob')
% xlabel('Correlator Output In-Phase Component');
% ylabel('Correlator Output Quadrature Component'), grid on
% legend('Noisy Observations','Model Prediction (Bayes)','Model Prediction (Myers)')
% % export_fig KF_est -eps -transparent

figure,
plot(1:N, (x-x_kf_est).^2);
xlabel('Sample');
ylabel('Squared State Error'), grid on
legend('Reciever Dynamic Phase Error','Carrier Doppler frequency shift Error','Carrier Doppler frequency rate Error')


figure,
subplot(2,1,1)
plot(Ntransient+1:N,R_RMSE(1,Ntransient+1:end),'r'), grid, hold on,
plot(Ntransient+1:N,R_RMSE3(1,Ntransient+1:end),'g'), grid,% hold on,
xlabel('Samples'), ylabel('RMSE($\sigma^2_{y_{k,i}}$)','Interpreter','latex')
legend('BCE','Myers'); %,'BCE (objective)','Myers')
title(sprintf('Covariance Estimation Error, Whiteness=[%f,%f]',whiteness(1),whiteness(2)))

subplot(2,1,2)
plot(Ntransient+1:N,R_RMSE(2,Ntransient+1:end),'r'), grid, hold on,
plot(Ntransient+1:N,R_RMSE3(2,Ntransient+1:end),'g'), grid,% hold on,
xlabel('Samples'), ylabel('RMSE($\sigma^2_{y_{k,q}}$)','Interpreter','latex')
legend('BCE', 'Myers'); %,'BCE (objective)','Myers')
% export_fig R_RMSE_case2 -eps -transparent

%     figure,
%     plot(Ntransient+1:N,mu_est(1,Ntransient+1:end), 'r'), grid, hold on,
%     plot(1:N,mean(yinnovations(1,:)).*ones(1,N), 'r--'),
%     plot(Ntransient+1:N,mu_est(2,Ntransient+1:end), 'b'),
%     plot(1:N,mean(yinnovations(2,:)).*ones(1,N), 'b-'), hold off
%     xlabel('Samples'), ylabel('$\bar{z}$','Interpreter','latex')
%     legend('Estimated Mean', 'Actual Mean');
%     title('Innovation Sequence Mean')
%



%plot

%end