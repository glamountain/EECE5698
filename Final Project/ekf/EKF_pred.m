function [x_kf_pred, P_kf_pred] = EKF_pred(f,Ff,Q,x_est_ant,P_est_ant)

% KF operating on a sample basis
% 
% v0: Pau Closas (06/2016)
% v1: Gerald LaMountain (04/2018)

Fk = Ff(x_est_ant);
x_kf_pred = f(x_est_ant);
P_kf_pred = Q + Fk*P_est_ant*Fk';



