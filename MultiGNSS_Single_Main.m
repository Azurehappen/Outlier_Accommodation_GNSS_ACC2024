% This code implement GNSS (GPS, GLO, GAL, BDS) under single frequency
% Special for CNES SSR data
% Correction type: PPP (Precise Point Positioning)
% clear all
% close all
%--------------------------------%
addpath('data')
addpath('parser')
addpath('time_compute')
addpath('eph')
addpath('pos')
addpath('corr')
%--------------------------------%
% Pick the Data Number
initpath = 'data/';
data_num = 1;
[eph_name,obs_name,IGS_name,data_base,code_bia,Grdpos,USTEC_folderpath] = datapathload(data_num,initpath);
%--------------------------------%

% Initialize parameters
[p,eph,obs] = initialization(eph_name,obs_name,Grdpos);
% obs.GLO.S1 = obs.GLO.S2; % Uncomment when GLO in data 8
% obs.BDS.P1(18,:) = 0;
% Mode setting
p.run_mode = 0;
p.post_mode  = 1; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
p.VRS_mode = 0;
p.IGS_enable = 1;
p.double_diff = 0;
p.elev_mark  = 15*pi/180;
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
p.inval = 1; % Computation time interval
p.tec_tmax = 15;
p.tec_tmin = 0;
%--------------------------------%
p = load_PPP_corr(p,data_base,IGS_name,eph,USTEC_folderpath);
if ~isempty(code_bia)
    p.code_bia = parser_bia(code_bia);
end
%-------------%
output = compute_gnss_ecef(p,eph,obs);
% p.min_sv = 6;
% output = compute_gnss_rcvr(p,eph,obs);

%%
figure
scatter(p.t,output.err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('ECEF positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

figure
scatter(p.t,output.hor_err,'.')
xtickformat('yyyy-MM-dd HH:mm:ss')
title('Horizontal positioning error')
xlabel('Local time')
ylabel('Error, unit: meter');grid on

total = output.sv_num_GPS + output.sv_num_GAL + output.sv_num_BDS;
figure
scatter(p.t,output.sv_num_GPS,'.')
hold on
scatter(p.t,output.sv_num_GAL,'.')
hold on
scatter(p.t,output.sv_num_BDS,'.')
hold on
scatter(p.t,total,'.')
title('total satellites been used')
legend('GPS','GAL','BDS','Total')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');grid on
% % 
% figure
% scatter(p.t,output.rover_clk/p.c,'.')
% title('Local bias')
% xlabel('Receiver time using GPS second');
% ylabel('Clock bias, seconds');grid on
% 

% figure
% subplot(311)
% scatter(p.t,output.ned_err(1,:),'.')
% title('North Error in NED');grid on;
% subplot(312)
% scatter(p.t,output.ned_err(2,:),'.')
% title('East Error in NED');grid on;
% subplot(313)
% scatter(p.t,output.ned_err(3,:),'.')
% title('Down Error in NED')
% xlabel('Receiver time using GPS second');grid on;
% figure
% for i=1:32
%     scatter(output.gpst,output.res_GPS(i,:),'.')
%     hold on
% end
% hold off
% grid on
% title('GPS residual')
% xlabel('Receiver time using GPS second');
% ylabel('Residual, unit: meter');
% figure
% scatter(output.gpst,output.sv_num_GPS,'.')
% title('Number of GPS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GLO,'.')
% title('Number of GLO satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_GAL,'.')
% title('Number of GAL satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.sv_num_BDS,'.')
% title('Number of BDS satellites')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on