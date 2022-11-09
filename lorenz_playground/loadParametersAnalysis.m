% This script saves the final parameters applied for all the other sripts
% into a single structure that is easy to read and modify in case of
% necessity 

%--------------------------------------------------------------------------
globalParameters = struct;

globalParameters.diameter_ganymede               = 70.373;
globalParameters.thickness_annulus               = 3;          % [pixels] thickness annulus for radial brightness plot
globalParameters.mx                              = 0.0246;     % field of view x dierction [arsec]
globalParameters.my                              = 0.0246;     % field of view y dierction [arsec]
globalParameters.A_mirror                        = 45238.9342; % [cm2] mirror HST
globalParameters.ganymede_assumed_brightness     = 1.3;        % [KReyleights]
globalParameters.RGanymedeCm                     = 2643.1e5;   % [cm] radius ganymede in

% giorno method specific parameters
globalParameters.number_of_bins_for_giornoMethod = 20;                             % number of sections to use for the application of giorno method
globalParameters.inner_radius_giorno   = globalParameters.diameter_ganymede/2;      % inner radius of the annulus as defined by giorno
globalParameters.outer_radius_giorno   = globalParameters.diameter_ganymede/2 *1.2; % outer radius of the annulus as defined by giorno

%--------------------------------------------------------------------------
% observation oe9z03010

observation3010.file = "./oe9z03010_flt.fits";
observation3010.name      = 'oe9z03010';
observation3010.x_pcenter  = 166;       % pixels
observation3010.y_pcenter  = 349;       % pixels
observation3010.poly_order = 3;
observation3010.best_n0    = 4600;      % 1/cm^2
observation3010.zc        = 4.95e-14 ;
observation3010.g_factor  = 7.3e-14;
observation3010.poly_order = 3;
observation3010.north_pole_angle_direction = 27.0; %deg defined anti-clockwise from horizonatal axis

% observation oe9z03011

observation3011.file = "./oe9z03011_flt.fits";
observation3011.name      = 'oe9z03011';
observation3011.x_pcenter  = 166;       % pixels
observation3011.y_pcenter  = 349;       % pixels
observation3011.poly_order = 3;
observation3011.best_n0    = 4300;      % 1/cm^2
observation3011.zc        = 4.95e-14 ;
observation3011.g_factor  = 7.3e-14;
observation3011.poly_order = 3;
observation3011.north_pole_angle_direction = 27.0; %deg defined anti-clockwise from horizonatal axis

% observation oe9z04011

observation4011.file = "./oe9z04011_flt.fits";
observation4011.name = 'oe9z04011';
observation4011.x_pcenter  = 158;       % pixels
observation4011.y_pcenter  = 361;       % pixels
%observation4011.poly_order = 2;
%observation4011.best_n0    = 5130;      % 1/cm^2
observation4011.poly_order = 3;
observation4011.best_n0    = 5700;      % 1/cm^2
observation4011.zc        = 4.95e-14 ;
observation4011.g_factor  = 7.3e-14;
observation4011.north_pole_angle_direction = 24.4; %deg defined anti-clockwise from horizonatal axis

% observation oe9z04010

observation4010.file = "./oe9z04010_flt.fits";
observation4010.name = 'oe9z04010';
observation4010.x_pcenter  = 158;       % pixels
observation4010.y_pcenter  = 361;       % pixels
%observation4010.poly_order = 2;
%observation4010.best_n0    = 5130;      % 1/cm^2
observation4010.poly_order = 3;
observation4010.best_n0    = 6800;      % 1/cm^2
observation4010.zc         = 4.95e-14 ;
observation4010.g_factor   = 7.3e-14;
observation4010.north_pole_angle_direction = 24.4; %deg defined anti-clockwise from horizonatal axis





