%% Figures
load('./../output_linear_NoGlo.mat')

%% 3D Error CDF Plot
N = length(output.GDOPLS);
fig1 = figure(1); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','p','^','h','>','o'];
legend_list=cell(output.LSSn,1);
for idx = 1:output.LSSn
    [f,x] = ecdf(output.err_LSS(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:65000; %1:200:N;
    legend_list{idx}=strcat('LSS \nu = ', num2str(output.LSSLambda(idx)));
end

leg1 = legend(legend_list); 
set(leg1,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of 3D Positioning Error');
xlabel('3D Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 2e0]); grid on

%% 3D Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 inf];
Marker = 'b.';

fig2 = figure(2); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.err_LSS(1,:),sz,output.err_LSS(1,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.err_LSS(2,:),sz,output.err_LSS(2,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

subplot(223); hold on; grid on
scatter(output.gpst,output.err_LSS(3,:),sz,output.err_LSS(3,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(224); hold on; grid on
scatter(output.gpst,output.err_LSS(4,:),sz,output.err_LSS(4,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

%% GDOP Scatter Plot
sz = 15; % size of the dots
GDOPLimits = [0 inf];
Marker = 'b.';

fig3 = figure(3); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.GDOPLSS(1,:),sz,output.GDOPLSS(1,:), Marker)
ylim(GDOPLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.GDOPLSS(2,:),sz,output.GDOPLSS(2,:), Marker)
ylim(GDOPLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(2))];
title(titleText,'Interpreter','tex')

subplot(223); hold on; grid on
scatter(output.gpst,output.GDOPLSS(3,:),sz,output.GDOPLSS(3,:), Marker)
ylim(GDOPLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.GDOPLSS(4,:),sz,output.GDOPLSS(4,:), Marker)
ylim(GDOPLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(4))];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')
sgtitle('GDOP')

%% Horizontal Error CDF Plot
N = length(output.GDOPLS);
fig4 = figure(4); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','>','p','h','^','o'];
legend_list=cell(output.LSSn,1);
for idx = 1:output.LSSn
    [f,x] = ecdf(output.hor_err_LSS(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('LSS \nu = ', num2str(output.LSSLambda(idx)));
end

leg4 = legend(legend_list); 
set(leg4,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of Horizontal Positioning Error');
xlabel('Horizontal Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 1e2]); grid on

%% Horizontal Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 1.5];
Marker = 'b.';

fig5 = figure(5); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.hor_err_LSS(1,:),sz,output.hor_err_LSS(1,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.hor_err_LSS(2,:),sz,output.hor_err_LSS(2,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(223); hold on; grid on
scatter(output.gpst,output.hor_err_LSS(3,:),sz,output.hor_err_LSS(3,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.hor_err_LSS(4,:),sz,output.hor_err_LSS(4,:), Marker)
ylim(errorLimits)
titleText = ['LSS, \nu = ',num2str(output.LSSLambda(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')

