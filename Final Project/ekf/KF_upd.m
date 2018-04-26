function [x_kf_est, P_kf_est, xinnovations] = KF_upd(yinnovations,H,R,x_kf_pred,P_kf_pred)

% KF operating on a sample basis
% 
% v0: Pau Closas (06/2016)

HPH = H*P_kf_pred*H';
S = HPH + R;
K_gain = P_kf_pred*H'/S;

x_kf_est = x_kf_pred + K_gain*yinnovations;
P_kf_est = P_kf_pred - K_gain*S*K_gain';

xinnovations = x_kf_est - x_kf_pred;