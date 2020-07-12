% $Date: $
% $Revision: $
% $Author: $
%=====================================================================
% This code implement GNSS (GPS, GLO, GAL, BDS) under single frequency
% Specially for Linear model
% Correction type: PPP (Precise Point Positioning)
clear all
% close all
%--------------------------------%
addpath('data')
addpath('parser')
addpath('time_compute')
addpath('eph')
addpath('pos')
addpath('corr')
addpath('linear')
addpath('linear/LIBRA')
%--------------------------------%
% Define the name of dataset from 'data' folder
% data 1
% dataname = 'COM3_190719_220647';
% % dataname = 'ppp_test';
% data_base = 'COM3_190719_220900';
% % data_base = [];

% data 2
% dataname = 'data/data2020/rover';
% data_base = 'data/data2020/base';
% % data_base = 'data/data2020/base_crfp';
% Grdpos = [-2430696.646;-4704190.763;3544329.081];

% data 3
dataname ='data/data202007/7-8-2020/COM3_200709_014237_run1';
data_base = 'data/data202007/7-8-2020/basedata';
M  = readmatrix('data/data202007/7-8-2020/grdtruth.csv');
Grdpos.pos = M(:,3:5)';
Grdpos.t = M(:,2);
%--------------------------------%
% Initialize parameters
[p,eph,obs] = initialization(dataname,dataname,Grdpos);
p.eph_b = eph;
% Measurement selection
p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
p.priorposcov = [2;2;2;1]; % cov([x0;clk;x_off])
p.ISBglo = 0.45; p.ISBglo_cov = 0.257^2; % Inter system bias GPS to GLO
p.ISBgal = 0.45; p.ISBgal_cov = 0.257^2; % Inter system bias GPS to GAL
p.ISBbds = 0.92; p.ISBbds_cov = 0.257^2; % Inter system bias GPS to BDS 
p.sig_y  = 0.72; % meters
% Mode setting
p.run_mode   = 1; %%%%% 0=real-time; 1=post processing
p.post_mode  = 2; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
p.elev_mask  = 10;
p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
p.enableGLO  = 0; % Enable GLO: 1 means enable, 0 means close
p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
p.enableBDS  = 0; % Enable BDS: 1 means enable, 0 means close
%--------------------------------%
% Enable Measurement selection
p.eb_LSS  = 1;
p.eb_MShb = 0;
p.eb_MStk = 0;
p.eb_LTS  = 1;
p.eb_TD   = 0;
%--------------------------------%
p.MShbConst = [0.25 0.5 1 1.345 1.75 2];
p.MStkConst = [2.5 3.0 3.5 4 4.5 4.685 5 5.5];
p.LTSOption = [1 2]; % 1 =default (check LTSlinear.m)
p.TDLambda  = [2.5 3 3.5 4 5];
p.LSSLambda = [0.8 1 1.2 1.5 2];
%--------------------------------%
p = load_PPP_corr(p,data_base,[],[],eph);
%--------------------------------%
% output = linear_gnss_ecef(p,eph,obs);
output = linear_gnss_ecef_parfor(p,eph,obs);
% save('output_linear_newlocal.mat','output','-v7.3')
%%
figure
scatter(output.gpst,output.hor_err_LS,'.')
title('LS ECEF error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');grid on

figure
scatter(output.gpst,output.hor_err_LSS(2,:),'.')
title('LSS ECEF error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');grid on

figure
scatter(output.gpst,output.hor_err_LTS(1,:),'.')
title('LTS ECEF error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');grid on

figure
scatter(output.gpst,output.err_MStk(6,:),'.')
title('MStk ECEF error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');grid on

figure
scatter(output.gpst,output.err_TD(3,:),'.')
title('TD ECEF error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');grid on

% figure
% scatter(output.gpst,output.hor_err_MShb,'.')
% title('Horizontal positioning error')
% xlabel('Receiver time using GPS second')
% ylabel('Error, unit: meter');grid on
% 
% figure
% scatter(output.gpst,output.hor_err_LTS,'.')
% title('Horizontal positioning error')
% xlabel('Receiver time using GPS second')
% ylabel('Error, unit: meter');grid on

% figure
% scatter(output.gpst,output.sv_num_GPS+output.sv_num_GLO+output.sv_num_GAL+output.sv_num_BDS,'.')
% title('total satellites been used')
% xlabel('Receiver time using GPS second')
% ylabel('Distance, unit: meter');grid on

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

% figure
% plot(output.gpst,output.cost,'.')
% title('Measurement residual comparison')
% xlabel('Receiver time using GPS second')
% title('Cost vs time')
% xlabel('Receiver time using GPS second')
% ylabel('Cost: $$\sum_{i=1}^{(m+n)} (\textbf{r}_i(\textbf{x}))^2$$',...
%     'Interpreter','Latex')