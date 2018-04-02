%
%   EECE5698-ST: GNSS signal processing
%       Simple pseudorange simulator to test navigation solvers
%
%   Pau Closas, Spring 2018 (using matlab functions from SW accompanying P. Groves book)

clearvars
close all
clc

rad2deg = 180/pi;
deg2rad = pi/180;


%% Simulation Parameters 

% number of epochs or pseudoranges
GNSS_config.no_epochs = 60;          
% time between consecutive pseudorange measurements (s)
GNSS_config.sampling = 0.01;          
% number of satellites
GNSS_config.no_sat = 30;              
% Mask angle (deg)
GNSS_config.mask_angle = 10;
% Code tracking error SD at zenith (m)
GNSS_config.code_track_err_SD = 2;
% Range rate tracking error SD at zenith (m/s)
GNSS_config.rate_track_err_SD = 0.02;
% Receiver clock offset at time=0 (m);
GNSS_config.rx_clock_offset = 10000;
% Receiver clock drift at time=0 (m/s);
GNSS_config.rx_clock_drift = 100;
% Initial estimated position (meters, ECEF)
GNSS_config.init_est_p_eb_ecef = [0;0;0];  

% specify location (static) - for instance: 42°20'13.4"N 71°05'25.9"W
true_phi_b_dms = [42 20 13.4];        % latitude (degree, minutes, seconds)
true_lambda_b_dms = -[71 05 25.9];    % longitude (degree, minutes, seconds)
true_h_b = 10;                        % height (meters)

% some transformations
true_phi_b_deg = true_phi_b_dms(1) + true_phi_b_dms(2)/60 + true_phi_b_dms(3)/3600;                % latitude (decimal degrees)
true_phi_b_rad = deg2rad*true_phi_b_deg;                                                                   % latitude (radians)

true_lambda_b_deg = true_lambda_b_dms(1) + true_lambda_b_dms(2)/60 + true_lambda_b_dms(3)/3600;    % longitude (decimal degrees)
true_lambda_b_rad = deg2rad*true_lambda_b_deg;                                                             % longitude (radians)

true_p_eb_ecef = pv_NED_to_ECEF(true_phi_b_rad,true_lambda_b_rad,true_h_b,[0;0;0]);                     % Ellipsoidal-to-Cartesian

% norm(true_p_eb_ecef)        % should be comparable to the Earh radius (R_0=6378 km < norm(position_ecef) < R_P=6356 km)

% % nice google map plot
figure, plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20)
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    plot_google_map

% user's velocity (Cartesian, ECEF)
true_v_eb_ecef = [0 0 0]';          % static receiver (m/s)

% initialize variables
time = 0;
dummy = 0;
pseudorange = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
pseudorangerate = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
elevations = zeros(GNSS_config.no_sat,GNSS_config.no_epochs);
est_p_eb_ecef_LS = zeros(3,GNSS_config.no_epochs);
est_clock_LS = zeros(1,GNSS_config.no_epochs);
est_phi_b_LS = zeros(1,GNSS_config.no_epochs);
est_lambda_b_LS = zeros(1,GNSS_config.no_epochs);
est_h_b_LS = zeros(1,GNSS_config.no_epochs);

%% main loop
for epoch = 1:GNSS_config.no_epochs    
    
    %% Generate pseudoranges

    % Determine satellite positions and velocities
    [sat_pos_es_e,sat_vel_es_e] = Satellite_positions_and_velocities(time,GNSS_config);
    % norm(sat_vel_es_e(1,:))        % magnitude of velocity for a satellite (m/s)

    % Generate GNSS measurements
    [GNSS_measurements,no_GNSS_meas] = Generate_GNSS_measurements(...
        time,sat_pos_es_e,sat_vel_es_e,true_p_eb_ecef,true_phi_b_rad,true_lambda_b_rad,true_v_eb_ecef,GNSS_config);

    pseudorange(:,epoch) = GNSS_measurements(:,1);          % for plotting
    pseudorangerate(:,epoch) = GNSS_measurements(:,2);
    elevations(:,epoch) = GNSS_measurements(:,9);
    
    %% Navigation solution
    % LS
    [est_p_eb_ecef_LS(:,epoch),est_clock_LS(epoch)] = GNSS_LS_position(GNSS_measurements,no_GNSS_meas,GNSS_config.init_est_p_eb_ecef);
    
    % transform estimates to latitude, longitude, and height
    [est_phi_b_LS(epoch),est_lambda_b_LS(epoch),est_h_b_LS(epoch),dummy] = pv_ECEF_to_NED(est_p_eb_ecef_LS(:,epoch),[0;0;0]);

    time = time + GNSS_config.sampling;
end

%% Compute error statistics (RMSE)
error_ecef_LS = sqrt(sum(mean((est_p_eb_ecef_LS - true_p_eb_ecef*ones(1,GNSS_config.no_epochs)).^2,2)))

%% Plot figures
t_vec = (0:GNSS_config.no_epochs-1)*GNSS_config.sampling;

% all pseudoranges
figure,
plot(t_vec,pseudorange(1:no_GNSS_meas,:).'), grid
    xlabel('time[s]'), ylabel('Pseudorange [m]')

% pick one satellite    
figure,
plot(t_vec,pseudorange(1,:).'), grid
    xlabel('time[s]'), ylabel('Pseudorange [m]')

% pseudorange(1,1) - pseudorange(1,end)    % pseudorange difference in no_epochs 
    
% elevations for all satellites    
figure,
plot(t_vec,rad2deg*elevations(1:no_GNSS_meas,:).'), grid
    xlabel('time[s]'), ylabel('elevations [m]')
    
% elevations for one satellites    
figure,
plot(t_vec,rad2deg*elevations(1,:).'), grid
    xlabel('time[s]'), ylabel('elevations [m]')

% true and estimated locations (LS)    
figure, plot(rad2deg*est_lambda_b_LS, rad2deg*est_phi_b_LS,'x','MarkerSize',5), hold on
    plot(true_lambda_b_deg, true_phi_b_deg,'.r','MarkerSize',20)
    xlabel('longitude [^o]'), ylabel('latitude [^o]')
    axis([true_lambda_b_deg-0.001 true_lambda_b_deg+0.001 true_phi_b_deg-0.001 true_phi_b_deg+0.001])
    plot_google_map    
    legend('LS estimates','True locations')
    
 