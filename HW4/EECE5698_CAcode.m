%   EECE5698-ST: GNSS signal processing
%       GPS L1 C/A code generation
%
%   Pau Closas, Spring 2018

clearvars
close all
clc

% satellite(s) ID
svnum = 19;
% C/A chip rate
Rc = 1.023e6;
% chip period
Tc = 1/Rc;

%% Generate and plot code(s)
[ca_code]=CAcodegen(svnum);
[ca_code2]=CAcodegen(svnum+1);

figure;
subplot(2,1,1); plot(ca_code);
subplot(2,1,2); plot(ca_code2);

%% Compute and plot autocorrelation(s)
ca_acor_line = xcorr (ca_code,ca_code);
ca_acor_circ = xcorr_circ(ca_code,ca_code);

figure;
subplot(2,1,1); plot(ca_acor_line);
subplot(2,1,2); plot(ca_acor_circ);


%% Compute and plot crosscorrelation with concatenated replicas
ca_ccor_line_cat = xcorr (ca_code,[ca_code ca_code ca_code]);
figure; plot(ca_ccor_line_cat);

%% Compute and plot crosscorrelation(s)
ca_xcor_line = xcorr (ca_code,ca_code2);
ca_xcor_circ = xcorr_circ(ca_code,ca_code2);

figure;
subplot(2,1,1); plot(ca_xcor_line);
subplot(2,1,2); plot(ca_xcor_circ);