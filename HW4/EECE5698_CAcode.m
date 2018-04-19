%   EECE5698-ST: GNSS signal processing
%       GPS L1 C/A code generation
%
%   Pau Closas, Gerald LaMountain, Spring 2018

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

fig1 = figure;
subplot(2,1,1); stem(ca_code);
title(sprintf('C/A Code for Satellite Vehicle #%d',svnum));
subplot(2,1,2); stem(ca_code2);
title(sprintf('C/A Code for Satellite Vehicle #%d',svnum+1));
saveas(fig1,['.\latex\figures\HW4P2a.png'])


%% Compute and plot autocorrelation(s)
[ca_acor_line, ll] = xcorr (ca_code,ca_code);
[ca_acor_circ, lc] = xcorr_circ(ca_code,ca_code);

fig2 = figure;
subplot(2,1,1); stem(ll, ca_acor_line);
title(sprintf('Linear Auto-Correlation of C/A Code for Satellite Vehicle #%d',svnum));
xlabel('Lag');
subplot(2,1,2); stem(lc, ca_acor_circ);
title(sprintf('Circular Auto-Correlation of C/A Code for Satellite Vehicle #%d',svnum));
xlabel('Lag');
saveas(fig2,['.\latex\figures\HW4P2b.png'])



%% Compute and plot crosscorrelation with concatenated replicas
[ca_ccor_line_cat, ll] = xcorr (ca_code,[ca_code ca_code ca_code]);
fig3 = figure; stem(ll, ca_ccor_line_cat);
title('Cross-Correlation between ca_code and thrice replicated c_code');
xlabel('Lag');
saveas(fig3,['.\latex\figures\HW4P2c.png'])


%% Compute and plot crosscorrelation(s)
[ca_xcor_line, ll] = xcorr (ca_code,ca_code2);
[ca_xcor_circ, lc] = xcorr_circ(ca_code,ca_code2);

fig4 = figure;
subplot(2,1,1); stem(ll, ca_xcor_line);
title(sprintf('Linear Cross-Correlation between C/A Code for Satellite Vehicle #%d and #%d',svnum,svnum+1));
subplot(2,1,2); stem(lc, ca_xcor_circ);
title(sprintf('Circular Cross-Correlation between C/A Code for Satellite Vehicle #%d and #%d',svnum,svnum+1));
saveas(fig4,['.\latex\figures\HW4P2d.png'])