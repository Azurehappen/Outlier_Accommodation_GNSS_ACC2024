%% Figures
load('./../output_linear_NoGlo.mat')

%% 3D Error CDF Plot
N = length(output.GDOPLS);
figure(1); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','^','p','h','>','o'];
legend_list=cell(output.TDn,1);
for idx = 1:output.TDn
    [f,x] = ecdf(output.err_TD(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('TD \lambda = ', num2str(output.TDLambda(idx)));
end

Legend = legend(legend_list); 
set(Legend,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of 3D Positioning Error');
xlabel('3D Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]); grid on

%% 3D Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 inf];
Marker = 'b.';

figure(2); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.err_TD(1,:),sz,output.err_TD(1,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.err_TD(2,:),sz,output.err_TD(2,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

subplot(223); hold on; grid on
scatter(output.gpst,output.err_TD(3,:),sz,output.err_TD(3,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(224); hold on; grid on
scatter(output.gpst,output.err_TD(4,:),sz,output.err_TD(4,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

%% GDOP Scatter Plot
sz = 15; % size of the dots
GDOPLimits = [0 inf];
Marker = 'b.';

figure(3); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.GDOPTD(1,:),sz,output.GDOPTD(1,:), Marker)
ylim(GDOPLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.GDOPTD(2,:),sz,output.GDOPTD(2,:), Marker)
ylim(GDOPLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(2))];
title(titleText,'Interpreter','tex')

subplot(223); hold on; grid on
scatter(output.gpst,output.GDOPTD(3,:),sz,output.GDOPTD(3,:), Marker)
ylim(GDOPLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.GDOPTD(4,:),sz,output.GDOPTD(4,:), Marker)
ylim(GDOPLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(4))];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')
sgtitle('GDOP')

%% Horizontal Error CDF Plot
N = length(output.GDOPLS);
figure(4); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','>','p','h','^','o'];
legend_list=cell(output.TDn,1);
for idx = 1:output.TDn
    [f,x] = ecdf(output.hor_err_TD(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('TD \lambda = ', num2str(output.TDLambda(idx)));
end

Legend = legend(legend_list); 
set(Legend,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of Horizontal Positioning Error');
xlabel('Horizontal Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 1e1]); grid on

%% Horizontal Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 1.5];
Marker = 'b.';

figure(5); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.hor_err_TD(1,:),sz,output.hor_err_TD(1,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.hor_err_TD(2,:),sz,output.hor_err_TD(2,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(223); hold on; grid on
scatter(output.gpst,output.hor_err_TD(3,:),sz,output.hor_err_TD(3,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.hor_err_TD(4,:),sz,output.hor_err_TD(4,:), Marker)
ylim(errorLimits)
titleText = ['Threshold Decisions, \lambda = ',num2str(output.TDLambda(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')
