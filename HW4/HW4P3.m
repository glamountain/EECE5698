clearvars
close all
clc

% Load hidden code
load('HW4P3.mat');
[~,len1] = size(ca_code_hidden);

% % Load all codes
% load('Cacodes.mat');

% Generate codes
ca_code_matrix = CAcodegen([1:32]);

[nc,len2] = size(ca_code_matrix);

xcor_mat = zeros(nc,len1+len2-1);
for i = 1:size(ca_code_matrix,1)
    xcor_mat(i,:) = xcorr(ca_code_hidden, ca_code_matrix(i,:));
end

% Find maximum
[~,sat] = max(max(xcor_mat,[],2));

% Plot cross correlations
fig1 = figure;
surf(xcor_mat)
title(sprintf('Hidden code cross correlation with known C/A Codes, Max=Sat#%d',sat))
xlabel('Lag'); ylabel('Sattelite Vehicle');
zlabel('Cross Correlation');
saveas(fig1,['.\latex\figures\HW4P3.png'])