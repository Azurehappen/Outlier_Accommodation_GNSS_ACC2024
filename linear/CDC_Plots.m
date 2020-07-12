% $Date: $
% $Revision: $
% $Author: $
%=====================================================================

% This code implement GNSS (GPS, GLO, GAL, BDS) under single frequency
% Specially for Linear model
% Correction type: PPP (Precise Point Positioning)
clear all
errorfname = 'linear\PlotsLinear\error1.png';
costfname = 'linear\PlotsLinear\cost1.png';
% close all
%--------------------------------%
addpath('./../data')
addpath('./../parser')
addpath('./../time_compute')
addpath('./../eph')
addpath('./../pos')
addpath('./../corr')
addpath('./../linear')
addpath('./../linear/LIBRA')
addpath('./../')
%%%=========================Least Squares==================================
%--------------------------------%
File = './output/output_CDC_old_struct.mat'; 
if ~isfile(File)
    disp('Now Performing Computations...')
    % Define the name of dataset from 'data' folder
    dataname = 'COM3_190719_220647';
    % dataname = 'ppp_test';
    data_base = 'COM3_190719_220900';
    % data_base = [];
    %--------------------------------%
    % Initialize parameters
    [p,eph,obs] = initialization(dataname);
    % Measurement selection
    p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
    p.priorposstd = [1.414;1.414;1.414;1]; % std of [x0;clk;x_off]
    p.ISBglo = 0.45; p.ISBglo_cov = 0.257; % Inner system bias GPS to GLO
    p.ISBgal = 0.45; p.ISBgal_cov = 0.257; % Inner system bias GPS to GAL
    p.ISBbds = 0.92; p.ISBbds_cov = 0.257; % Inner system bias GPS to BDS
    p.sig_y = 0.72; % meters
    p.select = 0; % 0=LS, 1=TD, 2=LSS, 3=MER, 4=LTS
    % Mode setting
    p.run_mode   = 1; %%%%% 0=real-time; 1=post processing
    p.post_mode  = 2; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
    p.elev_mask  = 10;
    p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
    p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
    p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
    p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
    %--------------------------------%
    p = load_PPP_corr(p,data_base,[],[]);
    
    %--------------------------------%
    LSoutput = linear_gnss_ecef(p,eph,obs);
    
    %%%======================Threshold Decisions===============================
    %--------------------------------%
    % Initialize parameters
    [p,eph,obs] = initialization(dataname);
    % Measurement selection
    p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
    p.priorposstd = [1.414;1.414;1.414;1]; % std of [x0;clk;x_off]
    p.ISBglo = 0.45; p.ISBglo_cov = 0.257; % Inner system bias GPS to GLO
    p.ISBgal = 0.45; p.ISBgal_cov = 0.257; % Inner system bias GPS to GAL
    p.ISBbds = 0.92; p.ISBbds_cov = 0.257; % Inner system bias GPS to BDS
    p.sig_y = 0.72; % meters
    p.select = 1; % 0=LS, 1=TD, 2=LSS, 3=MER, 4=LTS
    % Mode setting
    p.run_mode   = 1; %%%%% 0=real-time; 1=post processing
    p.post_mode  = 2; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
    p.elev_mask  = 10;
    p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
    p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
    p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
    p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
    %--------------------------------%
    p = load_PPP_corr(p,data_base,[],[]);
    
    %--------------------------------%
    TDoutput = linear_gnss_ecef(p,eph,obs);
    
    %%%=========================M-Estimator====================================
    %--------------------------------%
    % Initialize parameters
    [p,eph,obs] = initialization(dataname);
    % Measurement selection
    p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
    p.priorposstd = [1.414;1.414;1.414;1]; % std of [x0;clk;x_off]
    p.ISBglo = 0.45; p.ISBglo_cov = 0.257; % Inner system bias GPS to GLO
    p.ISBgal = 0.45; p.ISBgal_cov = 0.257; % Inner system bias GPS to GAL
    p.ISBbds = 0.92; p.ISBbds_cov = 0.257; % Inner system bias GPS to BDS
    p.sig_y = 0.72; % meters
    p.select = 3; % 0=LS, 1=TD, 2=LSS, 3=MER, 4=LTS
    % Mode setting
    p.run_mode   = 1; %%%%% 0=real-time; 1=post processing
    p.post_mode  = 2; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
    p.elev_mask  = 10;
    p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
    p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
    p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
    p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
    %--------------------------------%
    p = load_PPP_corr(p,data_base,[],[]);
    
    %--------------------------------%
    MERoutput = linear_gnss_ecef(p,eph,obs);
    
    %%%======================Least Trimmed Squares=============================
    %--------------------------------%
    % Initialize parameters
    [p,eph,obs] = initialization(dataname);
    % Measurement selection
    p.x0 = [-2430696.646;-4704190.763;3544329.081;0]; % prior state(1:3)
    p.priorposstd = [1.414;1.414;1.414;1]; % std of [x0;clk;x_off]
    p.ISBglo = 0.45; p.ISBglo_cov = 0.257; % Inner system bias GPS to GLO
    p.ISBgal = 0.45; p.ISBgal_cov = 0.257; % Inner system bias GPS to GAL
    p.ISBbds = 0.92; p.ISBbds_cov = 0.257; % Inner system bias GPS to BDS
    p.sig_y = 0.72; % meters
    p.select = 4; % 0=LS, 1=TD, 2=LSS, 3=MER, 4=LTS
    % Mode setting
    p.run_mode   = 1; %%%%% 0=real-time; 1=post processing
    p.post_mode  = 2; %%%% 0=Standard GNSS, 1 = PPP, 2= DGNSS
    p.elev_mask  = 10;
    p.enableGPS  = 1; % Enable GPS: 1 means enable, 0 means close
    p.enableGLO  = 1; % Enable GLO: 1 means enable, 0 means close
    p.enableGAL  = 1; % Enable GAL: 1 means enable, 0 means close
    p.enableBDS  = 1; % Enable BDS: 1 means enable, 0 means close
    %--------------------------------%
    p = load_PPP_corr(p,data_base,[],[]);
    
    %--------------------------------%
    LTSoutput = linear_gnss_ecef(p,eph,obs);
    
    %%%========================================================================
    
    % SAVE DATA SO WE DON'T HAVE TO RUN CODE EVERYTIME
    save(File, 'MERoutput', 'LTSoutput', 'TDoutput', 'LSoutput');
   
else 
    load(File);
end
%% Figures
%%%%===================Horizontal Positioning error========================
sz = 20; % size of the dots
costLimits = [0 250];
errorLimits = [0 4];
Marker = 'b.';

fig1 = figure(1); clf
colormap jet
subplot(221); hold on; grid on
scatter(LSoutput.gpst,LSoutput.hor_err,sz,LSoutput.hor_err, Marker)
ylim(errorLimits)
title('Least Squares ')
ylabel('Error, unit: meter')

subplot(223); hold on; grid on
scatter(LTSoutput.gpst,LTSoutput.hor_err,sz,LTSoutput.hor_err, Marker)
ylim(errorLimits)
title('Least Trimmed Squares')
ylabel('Error, unit: meter')
xlabel('Receiver time using GPS second')

subplot(222); hold on; grid on
scatter(TDoutput.gpst,TDoutput.hor_err,sz,TDoutput.hor_err, Marker)
ylim(errorLimits)
title('Threshold Decisions ')
ylabel('Error, unit: meter')

subplot(224); hold on; grid on
scatter(MERoutput.gpst,MERoutput.hor_err,sz,MERoutput.hor_err, Marker)
ylim(errorLimits)
title('M-Estimator')
ylabel('Error, unit: meter')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')

%saveas(fig1,errorfname);
%%%%=====================Cost vs Time======================================
fig2 = figure(2); clf
colormap jet
subplot(221); hold on; grid on
scatter(LSoutput.gpst,LSoutput.cost,sz,LSoutput.cost, Marker)
ylim(costLimits)
title('Least Squares')

subplot(223); hold on; grid on
scatter(LTSoutput.gpst,LTSoutput.cost,sz,LTSoutput.cost, Marker)
ylim(costLimits)
title('Least Trimmed Squares')
xlabel('Receiver time using GPS second')

subplot(222); hold on; grid on
scatter(TDoutput.gpst,TDoutput.cost,sz,TDoutput.cost, Marker)
ylim(costLimits)
title('Threshold Decisions')

subplot(224); hold on; grid on
scatter(MERoutput.gpst,MERoutput.cost,sz,MERoutput.cost, Marker)
ylim(costLimits)
title('M-Estimator')
xlabel('Receiver time using GPS second')
sgtitle('Cost vs time')

%saveas(fig2,costfname);

%% CDF PLOTS

%Cost

N = length(LSoutput.cost);
fig3 = figure(3); clf

[f,x] = ecdf(LSoutput.cost);
h_ls = plot(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 1;
h_ls.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(TDoutput.cost);
h_td = plot(x,f); hold on;
h_td.Color = 'g';
h_td.Marker = '>';
h_td.LineWidth = 1;
h_td.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(LTSoutput.cost);
h_lts = plot(x,f); hold on;
h_lts.Color = 'b';
h_lts.Marker = 'd';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(MERoutput.cost);
h_mer = plot(x,f); hold on;
h_mer.Color = 'k';
h_mer.Marker = 's';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

legend('LS', 'TD', 'LTS', 'MER');
title('')
sgtitle('Cumulative Distribution Of The Cost');
xlabel('Cost');
ylabel('Cumulative Probability');
xlim([0 100]);
grid on;

% Horizontal Error
fig4 = figure(4); clf

[f,x] = ecdf(LSoutput.hor_err);
h_ls = semilogx(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 1;
h_ls.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(TDoutput.hor_err);
h_td = semilogx(x,f); hold on;
h_td.Color = 'g';
h_td.Marker = '>';
h_td.LineWidth = 1;
h_td.MarkerIndices = 1:1e4:N;

[f,x] = ecdf(LTSoutput.hor_err);
h_lts = semilogx(x,f); hold on;
h_lts.Color = 'b';
h_lts.Marker = 'd';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(MERoutput.hor_err);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'k';
h_mer.Marker = 's';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;


legend('LS', 'TD', 'LTS', 'MER');
title('')
sgtitle('Cumulative Distribution Of The Horizontal Error');
xlabel('Horizontal Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 1e1]);
grid on;


