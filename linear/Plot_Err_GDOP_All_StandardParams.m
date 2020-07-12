%% Load processed DGNSS data
clear; clc;
load('./../output_linear_new2.mat')
% gpst = output.gpst/3600;

%%
% t = datetime(2020,5,24,17,00,00):seconds(30):datetime(2020,5,26,16,59,30);
% t = datetime(2020,4,21,22-7,13,24):seconds(1):datetime(2020,4,22,20-7,31,21);
% t(80183:80278)=[];
time = datetime(obs.tr_prime') - hours(7); % GNSS time -> datetime & UTC -> PST
figure(1); clf;
subplot(221); hold on; grid on
plot(time,obs.GPS.P1(2,:),'.')
xtickformat('HH:mm'); xtickangle(20)

subplot(222); hold on; grid on
plot(time,obs.GAL.P1(2,:),'.')
xtickformat('dd HH:mm')
xtickangle(30)

subplot(223); hold on; grid on
plot(time,obs.BDS.P1(2,:),'.')
xtickformat('HH:mm')
xtickangle(20)

subplot(224); hold on; grid on
plot(time,obs.GLO.P1(2,:),'.')
xtickformat('dd HH:mm')
xtickangle(10)

% xtickformat('yyyy-MM-dd HH:mm:ss')

% Choose index of parameters to create plots
%% indices of standard parameter values
LSSn_i  = 2;
MShbn_i = 3;
MStkn_i = 5;
TDn_i   = 2;
LTSn_i  = 2;

%% indices of best parameter values
LSSn_i  = 3;
MShbn_i = 3;
MStkn_i = 5;
TDn_i   = 3;
LTSn_i  = 2;


%% %%================== nsv, GDOP, error Scatter Plot ===========================
nsvMarker  = '.';   dnsvMarker = 's';
npriorMarker = '.'; dnpriorMarker = '.';
GDOPMarker = 'r.';
errMarker  = '.'; 
sz = 10; % size of the dots
% gpst = output.gpst;
% xlimits = [-inf inf]; %8.1*1e4];
nsv_ylimits  = [0 26];
nsv_marker_color = [0.5 0 0]; % brown
dnsv_marker_color = [0 0 0.5]; % brown
nprior_marker_color = [0.5 0 0]; % brown
dnprior_marker_color = [0 0.5 0]; % brown
err_marker_color = [0.55 0.55 0.55]; % gray
err_ylimits  = [0 5];
leftyaxis_color = 'k';
rightyaxis_color = [0.5 0 0]; % brown
gpst = datetime(output.tr_prime') - hours(7); % GNSS time -> datetime & UTC -> PST
xtickformat = 'mm-dd,HH'; tickangle = 0;
xlimits = [gpst(1) gpst(end)];

fig1 = figure(1); clf
sub1 = subplot(231); hold on; grid on
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvLS,sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvLS,sz,dnsvMarker);
set(sc4,'CData',dnsv_marker_color)
yyaxis left; ylim(err_ylimits)
sc2 = scatter(gpst,output.GDOPLS,sz,GDOPMarker); % set(sc2,'SizeData',10)
sc3 = scatter(gpst,output.err_LS,sz,errMarker);
set(sc3,'CData',err_marker_color)

sc5 = scatter(gpst,output.dnsvLS,sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
title('$\bf{LS}$','Interpreter','Latex')
datetick('x',xtickformat); xtickangle(tickangle)
xlim(xlimits)
ylabel('AGDOP & error (m)')

sub2 = subplot(232); hold on; grid on
txt = ['$\bf{TD}~(\lambda$ = ',num2str(output.TDLambda(TDn_i)),')'];
title(txt,'Interpreter','Latex')
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvTD(TDn_i,:),sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvTD(TDn_i,:),sz,dnsvMarker); 
set(sc4,'CData',dnsv_marker_color)
yyaxis left; ylim(err_ylimits)
sc2 = scatter(gpst,output.GDOPTD(TDn_i,:),sz,GDOPMarker);
sc3 = scatter(gpst,output.err_TD(TDn_i,:),sz,errMarker);
set(sc3,'CData',err_marker_color)

sc5 = scatter(gpst,output.dnpriorTD(TDn_i,:),sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
datetick('x',xtickformat); xtickangle(tickangle);
xlim(xlimits)

sub3 = subplot(233); hold on; grid on
txt = ['$\bf{LSS}~(\nu$ = ',num2str(output.LSSLambda(LSSn_i)),')'];
title(txt,'Interpreter','Latex')
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvLSS(LSSn_i,:),sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvLSS(LSSn_i,:),sz,dnsvMarker); 
set(sc4,'CData',dnsv_marker_color)
ylabel('No. of measurements')
yyaxis left; ylim(err_ylimits)
s2 = scatter(gpst,output.GDOPLSS(LSSn_i,:),sz,GDOPMarker);
sc3 = scatter(gpst,output.err_LSS(LSSn_i,:),sz,errMarker);
set(sc3,'CData',err_marker_color)

sc5 = scatter(gpst,output.dnpriorLSS(LSSn_i,:),sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
datetick('x',xtickformat); xtickangle(tickangle)
xlim(xlimits)


sub4 = subplot(234); hold on; grid on
title('$\bf{LTS}$','Interpreter','Latex')
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvLTS(LTSn_i,:),sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvLTS(LTSn_i,:),sz,dnsvMarker); 
set(sc4,'CData',dnsv_marker_color)
yyaxis left; ylim(err_ylimits)
sc2 = scatter(gpst,output.GDOPLTS(LTSn_i,:),sz,GDOPMarker);
sc3 = scatter(gpst,output.err_LTS(LTSn_i,:),sz,errMarker);
set(sc3,'CData',err_marker_color)

sc5 = scatter(gpst,output.dnpriorLTS(LTSn_i,:),sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
datetick('x',xtickformat); xtickangle(tickangle)
xlim(xlimits)
ylabel('AGDOP & error (m)')
xlabel('date, local time(hr)')

sub5 = subplot(235); hold on; grid on
txt = ['$\bf{MER-Tukey}$ (c = ',num2str(output.MStkConst(MStkn_i)),')'];
title(txt,'Interpreter','Latex')
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvMStk(MStkn_i,:),sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvMStk(MStkn_i,:),sz,dnsvMarker); 
set(sc4,'CData',dnsv_marker_color)
yyaxis left; ylim(err_ylimits)
sc2 = scatter(gpst,output.GDOPMStk(MStkn_i,:),sz,GDOPMarker);
sc3 = scatter(gpst,output.err_MStk(MStkn_i,:),sz,errMarker);
set(sc3,'CData',err_marker_color)
sc5 = scatter(gpst,output.dnpriorMStk(MStkn_i,:),sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
datetick('x',xtickformat); xtickangle(tickangle)
xlim(xlimits)
xlabel('date, local time(hr)')

sub6 = subplot(236); hold on; grid on
txt = ['$\bf{MER-Huber}$ (c = ',num2str(output.MShbConst(MShbn_i)),')'];
title(txt,'Interpreter','Latex')
yyaxis right; ylim(nsv_ylimits)
sc1 = scatter(gpst,output.nsvMShb(MShbn_i,:),sz,nsvMarker); 
set(sc1,'CData',nsv_marker_color)
sc4 = scatter(gpst,output.dnsvMShb(MShbn_i,:),sz,dnsvMarker); 
set(sc4,'CData',dnsv_marker_color)
ylabel('No. of measurements')
yyaxis left; ylim(err_ylimits)
sc2 = scatter(gpst,output.GDOPMShb(MShbn_i,:),sz,GDOPMarker);
sc3 = scatter(gpst,output.err_MShb(MShbn_i,:),sz,errMarker);
set(sc3,'CData',err_marker_color)
sc5 = scatter(gpst,output.dnpriorMShb(MShbn_i,:),sz,dnpriorMarker); 
set(sc5,'CData',dnprior_marker_color)
ax = gca; 
ax.YAxis(1).Color = leftyaxis_color; 
ax.YAxis(2).Color = rightyaxis_color;
datetick('x',xtickformat); xtickangle(tickangle)
xlim(xlimits)
xlabel('date, local time(hr)')


Legend6 = legend({['Geometric Dilution of Precision (m)'  ''],...
    ['ECEF Positioning error (m)'],...
    ['No. of unused measurements'], ['No. of measurements used']},'FontSize',9);

sub1.Position(1) =  sub1.Position(1) - 0.08;
sub2.Position(1) =  sub2.Position(1) - 0.1;
sub3.Position(1) =  sub3.Position(1) - 0.12;
sub4.Position(1) =  sub4.Position(1) - 0.08; sub4.Position(2) = sub4.Position(2) + 0.02;
sub5.Position(1) =  sub5.Position(1) - 0.1;  sub5.Position(2) = sub5.Position(2) + 0.02;
sub6.Position(1) =  sub6.Position(1) - 0.12; sub6.Position(2) = sub6.Position(2) + 0.02;
set(gcf,'units','inches','position',[0,2,16,6]) % position of figure window
set(Legend6,'Orientation','Horizontal');
Legend6.Position = [0.200,0.020 - 0.01,0.427,0.033];

%% Posterior Standard Deviation
N = size(output.err_LS(1,:),2)
std_LS   = nanstd(output.err_LS(1,:))
std_LSS  = nanstd(output.err_LSS(LSSn_i,:))
std_LTS  = nanstd(output.err_LTS(1,:))
std_MStk = nanstd(output.err_MStk(MStkn_i,:))
std_MShb = nanstd(output.err_MShb(MShbn_i,:))
std_TD   = nanstd(output.err_TD(TDn_i,:))

%% ======================= 3D Error CDF only ==============================
N = length(output.GDOPLS);
MarkerIndices = 1:5e3:65000;
% colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
% markerlist = ['o','h','d','v','s','p','x','<','^','h','>','o'];

fig2 = figure(2); clf
[f,x] = ecdf(output.err_LS);
h_ls = semilogx(x,f); hold on;
h_ls.Color = 'r';
h_ls.Marker = 'o';
h_ls.LineWidth = 4;
h_ls.LineWidth = 1;
h_ls.MarkerIndices = MarkerIndices;

[f,x] = ecdf(output.err_LSS(LSSn_i,:));
h_lss = semilogx(x,f); hold on;
h_lss.Color = 'k';
h_lss.Marker = 'p';  
h_lss.MarkerSize = 3;
h_lss.LineWidth = 1; 
h_lss.MarkerIndices = MarkerIndices;

[f,x] = ecdf(output.err_LTS);
h_lts = semilogx(x,f); hold on;
h_lts.Color = 'g';
h_lts.Marker = '*';
h_lts.LineWidth = 1;
h_lts.MarkerIndices = MarkerIndices;

[f,x] = ecdf(output.err_MShb(MShbn_i,:));
h_mshb = semilogx(x,f); hold on;
h_mshb.Color = 'b';
h_mshb.Marker = 'd';
h_mshb.MarkerSize = 4;
h_mshb.LineWidth = 1;
h_mshb.MarkerIndices = MarkerIndices;

[f,x] = ecdf(output.err_MStk(MStkn_i,:));
h_mstk = semilogx(x,f); hold on;
h_mstk.Color = 'm';
h_mstk.Marker = 's';
h_mstk.MarkerSize = 7;
h_mstk.LineWidth = 1;
h_mstk.MarkerIndices = MarkerIndices;

[f,x] = ecdf(output.err_TD(TDn_i,:));
h_td = semilogx(x,f); hold on;
h_td.Color = 'c';
h_td.Marker = '^';
h_td.MarkerSize = 9;
h_td.LineWidth = 1;
h_td.MarkerIndices = MarkerIndices;

leg2 = legend('LS (no outlier removed)',['LSS (\nu = ',...
    num2str(output.LSSLambda(LSSn_i)),')'],'LTS',['MER-Huber (c = ',...
    num2str(output.MShbConst(MShbn_i)),')'],['MER-Tukey (c = ',...
    num2str(output.MStkConst(MStkn_i)),')'],['TD (\lambda = ',...
    num2str(output.TDLambda(TDn_i)),')']);
set(leg2,'Location','northwest');
set(leg2,'Interpreter','tex');
%sgtitle('Cumulative Distribution of 3D Positioning Error');
xlabel('Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 1e1]);
grid on;

%% %%=================== GDOP only Scatter Plot =========================

sz = 15; % size of the dots
Marker = 'b.';
GDOPaxislimits = [-inf inf -inf 3.2];

fig2 = figure(2); clf
%sgtitle('Geometric Dilution of Precision')
colormap jet
subplot(231); hold on; grid on
scatter(gpst,output.GDOPLS,sz,output.GDOPLS, Marker)
ylabel('GDOP, unit: meters');
title('$\bf{LS}$','Interpreter','Latex')
axis(GDOPaxislimits)

subplot(232); hold on; grid on
txt = ['$\bf{TD}~(\lambda$ = ',num2str(output.TDLambda(TDn_i)),')'];
title(txt,'Interpreter','Latex')
scatter(gpst,output.GDOPTD(TDn_i,:),sz,output.GDOPTD(TDn_i,:), Marker)
axis(GDOPaxislimits)

subplot(233); hold on; grid on
txt = ['$\bf{LSS}~(\nu$ = ',num2str(output.LSSLambda(LSSn_i)),')'];
title(txt,'Interpreter','Latex')
scatter(gpst,output.GDOPLSS(LSSn_i,:),sz,output.GDOPLSS(LSSn_i,:), Marker)
axis(GDOPaxislimits)

subplot(234); hold on; grid on
title('$\bf{LTS}$','Interpreter','Latex')
scatter(gpst,output.GDOPLTS,sz,output.GDOPLTS, Marker)
axis(GDOPaxislimits)
ylabel('GDOP, unit: meters'); xlabel('GPS second')

subplot(235); hold on; grid on
txt = ['$\bf{MER-Tukey}$ (c = ',num2str(output.MStkConst(MStkn_i)),')'];
title(txt,'Interpreter','Latex')
scatter(gpst,output.GDOPMStk(MStkn_i,:),sz,output.GDOPMStk(MStkn_i,:), Marker)
axis(GDOPaxislimits)
xlabel('GPS second')

subplot(236); hold on; grid on
txt = ['$\bf{MER-Huber}$ (c = ',num2str(output.MShbConst(MShbn_i)),')'];
title(txt,'Interpreter','Latex')
scatter(gpst,output.GDOPMShb(MShbn_i,:),sz,output.GDOPMShb(MShbn_i,:), Marker)
axis(GDOPaxislimits)
xlabel('GPS second')

set(gcf,'units','inches','position',[0,3,16,5])

%% %% 3D Positioning error only Scatter Plot
sz = 15; % size of the dots
axeslimits = [-inf inf 0 inf];
errMarker = '.';

fig3 = figure(3); clf
sgtitle('3D positioning error')
colormap jet
subplot(231); hold on; grid on
scatter(gpst,output.err_LS,sz,output.err_LS, errMarker)
axis(axeslimits)
titleText = ['Least Squares'];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(232); hold on; grid on
scatter(gpst,output.err_LSS(LSSn_i,:),sz,output.err_LSS(LSSn_i,:), errMarker)
axis(axeslimits)
titleText = ['LSS, \lambda = ',num2str(output.LSSLambda(LSSn_i))];
title(titleText,'Interpreter','tex')

subplot(233); hold on; grid on
scatter(gpst,output.err_TD(TDn_i,:),sz,output.err_TD(TDn_i,:), errMarker)
axis(axeslimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(TDn_i))];
title(titleText,'Interpreter','tex')

subplot(234); hold on; grid on
scatter(gpst,output.err_MShb(MShbn_i,:),sz,output.err_MShb(MShbn_i,:), errMarker)
axis(axeslimits)
titleText = ['M-Estimator(Huber), Const = ',num2str(output.MShbConst(MShbn_i))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')

subplot(235); hold on; grid on
scatter(gpst,output.err_MStk(MStkn_i,:),sz,output.err_MStk(MStkn_i,:), errMarker)
axis(axeslimits)
titleText = ['M-Estimator (Tukey), BisqConst = ',num2str(output.MStkConst(MStkn_i))];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')

subplot(236); hold on; grid on
scatter(gpst,output.err_LTS(1,:),sz,output.err_LTS(1,:), errMarker)
axis(axeslimits)
titleText = ['LTS'];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')


