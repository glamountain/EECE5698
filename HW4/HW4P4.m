clearvars
close all
clc

% signal parameters
fs = 10e6;          % Sampling frequency
fIF = 0;            % Intermediate frequency
fc = 1.023e6;       % Code rate of the code [Hz]
Tcoh = 0.001;       % Coherent integration time [s]
Nc = Tcoh * fs;     % number of samples contained in the coherent integration time

% Acquisition parameters
Nd = 81;            % Number of Doppler bins
DopStep = 125;      % Doppler bin size in Hz
secondOfData = 0.1; % Seconds of data to read
Nnci = 10;           % number of non-coherent integrations of the CAF
SVIDs = 1:32;          % satellite vehicles (SV) to detect [7 16 19 21 22 25]

% read file with the IF capture
fid = fopen ('realGPSL1capture.bin','r');
[data, cnt_data] = fread(fid, 2 * secondOfData * fs, 'int8');
data = data(1:2:end) + 1i * data(2:2:end);

CAF_aux = zeros(Nd,Nc);
DopplerEst = -ones(1,32);
DelayEst = -ones(1,32);

% generate local replica and resample from Tc to Ts>Tc

% loop over all possible satellites
for svnum = SVIDs

    [ca_code]=CAcodegen(svnum);
    % Resample the code at data sampling frequency
    ca_code_resampled = ResampleCode( ca_code, Nc, fs, 0, fc );   
    
    CAF = 0;        % initialized CAF to zero every time
    
    % loop to average over noncoherent integrations
    for ii = 1:Nnci
        y =  data( (ii - 1) * Nc + (1:Nc) ).';   % use just 1 period of code at the time
        
        % loop over frequency bins
        for ff=1:Nd
            fdbin = fIF/fs + (ff - ceil(Nd/2))*DopStep/fs;   % normalized Doppler bin
            CAF_aux(ff,:) = computeCAF(y, ca_code_resampled, fdbin);
        end
        % integrate non-coherently the ii test statistics |CAF|^2 
        CAF = CAF + abs(CAF_aux).^2;

        % plot 2D grid search
        CAF_normalized = CAF/max(max(CAF));     % normalize to 1 the maximum value
        
    end
    
        figure(svnum)
        surf((0:(Nc - 1)) / fs, ((1:Nd) - ceil(Nd/2))*DopStep, CAF_normalized, 'EdgeColor', 'none');
        axis tight, set( gca, 'FontSize', 16 )
        xlabel('Code delay [s]'), ylabel('Doppler [Hz]')
        title(['SV ' num2str(svnum) ' , ' num2str(ii) ' non-coherent integrations'])
%         pause
    % pause
    
    % estimate Doppler (if satellite detected)
    [~, DopInd] = max(max(CAF.'));
    DopplerEst(svnum) = fIF + (DopInd - ceil(Nd/2))*DopStep;

    % estimate time-delay (if satellite detected)
    [~, codInd] = max(max(CAF));
    DelayEst(svnum) =  (codInd - 1) / fs;

end
