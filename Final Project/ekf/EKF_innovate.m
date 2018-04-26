function [yinnovations, HPH] = EKF_innovate(y,h,Hh,x_kf_pred,P_kf_pred)

% KF operating on a sample basis
% 
% v0: Gerald LaMountain (04/2018)

Hk  = Hh(x_kf_pred);
HPH = Hk*P_kf_pred*Hk';
yinnovations = y - h(x_kf_pred);



