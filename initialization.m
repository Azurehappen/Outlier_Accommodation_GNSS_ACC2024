function [p,eph,obs,t] = initialization(eph_name,obs_name,Grdpos)
%-------------------------------------------------------------------------%
%Define the constant parameters for GNSS system
%-------------------------------------------------------------------------%
% dataname: the rover data file name
%
%
%
%
%----------------------%
p.gps.num_prn = 33; % The amoount of GPS satellites
p.gal.num_prn = 38; % The amoount of GAL satellites
p.glo.num_prn = 34; % The amoount of GLO satellites
p.bds.num_prn = 50; % The amoount of BDS satellites
%----------------------%

matname = [eph_name '_nav.mat'];
navname = [eph_name '.nav'];
obsname = [obs_name '.obs'];
if exist(matname,'file')==2 % Check if the data already been parsed
    load(matname);
else
    % Get ephemeris data (.nav file, RINEX verion 3.03)
    eph = parser_nav(p,navname);
    save ([eph_name,'_nav.mat'], 'eph');
end
matname = [obs_name '_obs.mat'];
if exist(matname,'file')==2 % Check if the data already been parsed
    load(matname);
else
    % Get observables data (.obs file, RINEX verion 3.03)
    obs = parser_obs(obsname);
    save ([obs_name,'_obs.mat'], 'obs');
end
p.t = datetime(obs.tr_prime');
% load('data/DCB_GLO.mat');
% p.icb_glo = DCB_P1C1;
% Setting
p.inval = 1; % Computation time interval
p.run_mode=1; %%%%% 0=real-time; 1=post processing
p.post_mode = 1;%%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
p.iono_map=0; %%%%% 0=USTEC; 1=IGS ionex
p.ins=1;    %%%%%%% 0=no imu data; 1=use simulated/real imu data
p.imu_freq=200; %%%%%% frequency rate for simulated/real IMU data
p.P_base = [-2430697.667;-4704189.148;3544329.063]; % Base station position
p.Grdpos = Grdpos; % Ground truth of rover
p.freq = 1; %%%% 1 = single frequency, 2 = dual frequency
p.Ek0 = 0; % Initial condition of Ek
p.enableGPS = 0; % Enable GPS: 1 means enable, 0 means close
p.enableGLO = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL = 0; % Enable GAL: 1 means enable, 0 means close
p.enableBDS = 0; % Enable BDS: 1 means enable, 0 means close
p.IGS_enable = 1; % Enable IGS correction: 1 means enable, 0 means close
p.L1freq = 1575.42e6; % L1 frequency (Hz)
p.L2freq = 1227.6e6; % L2 frequency (Hz)
p.L1glo = 1602.0e6; % GLO L1 frequency (Hz)
p.L2glo = 1246.0e6; % GLO L2 frequency (Hz)
p.E1freq = 1575.42e6; % E1 frequency (Hz)
p.E6freq = 1278.75e6; % E6 frequency (Hz)
p.E5freq = 1191.795e6; % E5 frequency (Hz)
p.E5afreq = 1176.45e6; % E5a frequency (Hz)
p.E5bfreq = 1207.14e6; % E5b frequency (Hz)
p.B1freq = 1561.098e6; % B1I frequency (Hz)
p.B2afreq = 1207.14e6; % B2b frequency (currently not supported)
p.state0 = [0;0;0;0]; % Initial state vector at first iteration, maybe changed in other functions
%                             [x;y;z;clk_bias]
p.mk=0;
p.lat = 0; p.lon = 0; p.h_r = 0; % Initialize latitude, longitude and height
p.USTEC = 0;
p.eph_base = [];
p.obs_base = [];
p.IGS = [];
% Common parameters
p.c = 2.99792458e8; % Speed of light (m/s)
p.pi = pi; % Archimedes' constant
if eph.LeapSeconds~=0
    p.gpstime_ahead = eph.LeapSeconds; % The difference between GPS time and UTC, GPS has 18 seconds ahead of UTC
else
    % In some case, there is no leap seconds in eph data.
    % This constant need to be manually changed if Leap seconds change.
    p.gpstime_ahead = 18;
end
p.Re        = 6378136.3; %%%% radius of earth in m
p.h_iono    = 350000;    %%%% height of ionosphere of maximum TEC
% Earth model
p.g = 9.80665;                    % standard gravity
p.omge = 7.2921151467e-5;     % earth angular velocity (rad/s)
p.a = 6378137.0;                  % semi-major axis length, meters
p.b = 6356752.31424518;           % semi-minor axis length, meters
p.f = 1.0/298.257223563;          % flatness of ellipsoid
p.e = sqrt((p.a^2-p.b^2)/(p.a^2));      % first eccentricity of ellipsoid
p.ep = sqrt((p.a^2-p.b^2)/(p.b^2));     % second eccentricity of ellipsoid
%-------------------------------------------------------------------------%
% Treshold
p.elev_mark = 10*pi/180; % Elevation treshold
p.sig_strg = 0; % Signal strength treshold
p.satdelta = 1e-6; % In function 'eph2pos', the threshold for sat pos convergence
p.min_sv = 5; % the minimum of the number of satellites
p.LSthrsh = 1e-8; % threshold of delta_x in LS solver
% Iteration Number
p.NsatEk = 20; % In function 'eph2pos', the maxinum of iteration for Ek convergence
p.Nls = 20; % Least square
p.GDOP_mark = 30;
%-------------------------------------------------------------------------%
% Time Synchronization
% GPS starting at 1980-1-6 00:00:00
% Galileo starting at 1999-8-22 00:00:13
% BDS starting at 2006-1-1 00:00:00
% Leap seconds between 1980 and 1999 is 13s, Hence GPS seconds = GAL seconds
p.gal.lps_gps = 0;
% Leap seconds between 1980 and 2060 is 14s, Hence GPS seconds = BDS seconds + 14
p.bds.lps_gps = 14;
% GPS has 18 seconds ahead of UTC, GLONASS is synchronized to UCT
% Hence GPS seconds = GLO seconds + p.gpstime_ahead
p.glo.lps_gps = p.gpstime_ahead;

%---------------------------------------%  
% GPS Constants
% Reference 'https://www.gps.gov/technical/icwg/IS-GPS-200H.pdf'
%---------------------------------------%  
p.gps.mu = 3.986005e+14; % WGS84 value of the earth's gravitational constant for GPS user (m^3/s^2)
p.gps.OmegaDot_e = 7.2921151467e-5; % Earth rotation rate (rad/s)
p.gps.F = -4.442807633e-10; % IS-GPS-200H page 96
p.gps.message_duration = 7200; % Maximum time difference between obs data and eph message
% Galileo Constants
% Reference 'https://www.gsc-europa.eu/sites/default/files/sites/all/files/Galileo-OS-SIS-ICD.pdf'
%---------------------------------------%
p.gal.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.gal.OmegaDot_e = 7.2921151467e-5; % mean angular velocity of the Earth
p.gal.F = -4.442807309e-10; % Galileo-OS-SIS-ICD page 58
p.gal.message_duration = 7200;
% GLONASS Constants
% Reference ''
%---------------------------------------%
p.glo.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.glo.OmegaDot_e = 7.2921151467e-5; % Earth's rotation rate
p.glo.F = -2*sqrt(p.gal.mu)/(p.c^2); % 
p.glo.a_e = 6378136; %  semi-major (equatorial) axis of the PZ-90 Earth’s ellipsoid
p.glo.C_20 = 1082625.75e-9; %  second degree zonal coefficient of normal potential
p.glo.message_duration = 1800; %54000;
% BeiDou Constants
% Reference 'http://en.beidou.gov.cn/SYSTEMS/ICD/201902/P020190227702348791891.pdf'
%---------------------------------------%
p.bds.mu = 3.986004418e+14; % Geocentric gravitational constant (m^3/s^2)
p.bds.OmegaDot_e = 7.2921150e-5; % Earth's rotation rate
p.bds.F = -2*sqrt(p.gal.mu)/(p.c^2); % BeiDou-ICD page 58
p.bds.num_prn = 40; % The amoount of BDS satellites
p.bds.message_duration = 3600;
% Measurement selection
p.select = 0;

%---------------------------------------%
p.GPS_C1C = 1;p.GPS_C1W = 2;p.GPS_C2L = 3;p.GPS_C2W = 4;
p.GLO_C1C = 1;p.GLO_C1P = 2;p.GLO_C2C = 3;p.GLO_C2P = 4;
p.GAL_C1X = 1;p.GAL_C7X = 2;
P.BDS_C2I = 1;p.BDS_C7I = 2;

end
