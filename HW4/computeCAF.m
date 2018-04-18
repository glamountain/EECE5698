function CAF = computeCAF(y, ca_code, fdbin)

% Compute the CAF with input signal y[n], local code ca_code, at the Doppler frequency
% bin fdbin
%
%   EECE5698-ST: GNSS signal processing
%       GPS L1 C/A code generation
%
%   Pau Closas, Spring 2018

Nc = length(y);  % number of samples coherently integrated
n = 0:(Nc - 1);  % time index

y_freq = fft(y.*exp(-1i*2*pi*fdbin.*n));
ca_code_freq = fft(ca_code);
CAF = ifft(y_freq.*conj(ca_code_freq)) / Nc;