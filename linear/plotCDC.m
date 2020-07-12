%% Figures
load('output_CDC.mat')
N = length(output.gpst);
%%
%%%%===================Horizontal Positioning error========================
sz = 20; % size of the dots
costLimits = [0 250];
errorLimits = [0 4];
Marker = 'b.';

fig1 = figure(1); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.hor_err_LS,sz,output.hor_err_LS, Marker)
ylim(errorLimits)
title('Least Squares ')
ylabel('Error, unit: meter')

subplot(222); hold on; grid on
scatter(output.gpst,output.hor_err_LTS,sz,output.hor_err_LTS, Marker)
ylim(errorLimits)
title('Least Trimmed Squares')
ylabel('Error, unit: meter')
xlabel('Receiver time using GPS second')

subplot(223); hold on; grid on
scatter(output.gpst,output.hor_err_MShb,sz,output.hor_err_MShb, Marker)
ylim(errorLimits)
title('M-Estimator (Huber)')
ylabel('Error, unit: meter')

subplot(224); hold on; grid on
scatter(output.gpst,output.hor_err_MStk,sz,output.hor_err_MStk, Marker)
ylim(errorLimits)
title('M-Estimator (Tukey)')
ylabel('Error, unit: meter')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')

%%
%%%%=====================Cost vs Time======================================
% fig2 = figure(2); clf
% colormap jet
% subplot(221); hold on; grid on
% scatter(LSoutput.gpst,LSoutput.cost,sz,LSoutput.cost, Marker)
% ylim(costLimits)
% title('Least Squares')
% 
% subplot(223); hold on; grid on
% scatter(LTSoutput.gpst,LTSoutput.cost,sz,LTSoutput.cost, Marker)
% ylim(costLimits)
% title('Least Trimmed Squares')
% xlabel('Receiver time using GPS second')
% 
% subplot(222); hold on; grid on
% scatter(TDoutput.gpst,TDoutput.cost,sz,TDoutput.cost, Marker)
% ylim(costLimits)
% title('Threshold Decisions')
% 
% subplot(224); hold on; grid on
% scatter(MERoutput.gpst,MERoutput.cost,sz,MERoutput.cost, Marker)
% ylim(costLimits)
% title('M-Estimator')
% xlabel('Receiver time using GPS second')
% sgtitle('Cost vs time')

%saveas(fig2,costfname);

%% CDF PLOTS
% Cost CDF
fig3 = figure(3); clf

[f,x] = ecdf(output.costLS);
h_ls = semilogx(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 1;
h_ls.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costLTS);
h_td = semilogx(x,f); hold on;
h_td.Color = 'g';
h_td.Marker = '>';
h_td.LineWidth = 1;
h_td.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costMShb);
h_lts = semilogx(x,f); hold on;
h_lts.Color = 'b';
h_lts.Marker = 'd';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costMStk);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'k';
h_mer.Marker = 's';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costTD1);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'm';
h_mer.Marker = 'x';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costTD2);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'c';
h_mer.Marker = 'p';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.costTD3);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'y';
h_mer.Marker = 'h';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

legend('LS', 'LTS', 'MER Huber', 'MER Tukey', 'TD, L=1.0','TD, L=2.0','TD, L=3.0');
title('')
sgtitle('Cumulative Distribution Of The Cost');
xlabel('Cost');
ylabel('Cumulative Probability');
xlim([0.01 100]);
grid on;

%% Horizontal Error
fig4 = figure(4); clf

[f,x] = ecdf(output.hor_err_LS);
h_ls = semilogx(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 1;
h_ls.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_LTS);
h_td = semilogx(x,f); hold on;
h_td.Color = 'g';
h_td.Marker = '>';
h_td.LineWidth = 1;
h_td.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_MShb);
h_lts = semilogx(x,f); hold on;
h_lts.Color = 'b';
h_lts.Marker = 'd';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_MStk);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'k';
h_mer.Marker = 's';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_TD1);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'm';
h_mer.Marker = 'x';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_TD2);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'c';
h_mer.Marker = 'p';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.hor_err_TD3);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'y';
h_mer.Marker = 'h';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

legend('LS', 'LTS', 'MER Huber', 'MER Tukey', 'TD, L=1.0','TD, L=2.0','TD, L=3.0');
title('')
sgtitle('Cumulative Distribution Of The Horizontal Error');
xlabel('Horizontal Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]);
grid on;

%% 3D ECEF Error
fig4 = figure(4); clf

[f,x] = ecdf(output.err_LS);
h_ls = semilogx(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 1;
h_ls.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_LTS);
h_td = semilogx(x,f); hold on;
h_td.Color = 'g';
h_td.Marker = '>';
h_td.LineWidth = 1;
h_td.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_MShb);
h_lts = semilogx(x,f); hold on;
h_lts.Color = 'b';
h_lts.Marker = 'd';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_MStk);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'k';
h_mer.Marker = 's';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_TD1);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'm';
h_mer.Marker = 'x';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_TD2);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'c';
h_mer.Marker = 'p';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

[f,x] = ecdf(output.err_TD3);
h_mer = semilogx(x,f); hold on;
h_mer.Color = 'y';
h_mer.Marker = 'h';
h_mer.LineWidth = 1;
h_mer.MarkerIndices = 1:5e3:N;

legend('LS', 'LTS', 'MER Huber', 'MER Tukey', 'TD, L=1.0','TD, L=2.0','TD, L=3.0');
title('')
sgtitle('Cumulative Distribution Of The 3D Error');
xlabel('3D Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]);
grid on;