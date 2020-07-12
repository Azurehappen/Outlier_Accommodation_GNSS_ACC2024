%% Figures
load('./../output_linear.mat')

%% 3D Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 4];
Marker = 'b.';

fig1 = figure(1); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.err_LTS(1,:),sz,output.err_LTS(1,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.err_LTS(2,:),sz,output.err_LTS(2,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

subplot(223); hold on; grid on
scatter(output.gpst,output.err_LTS(3,:),sz,output.err_LTS(3,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(224); hold on; grid on
scatter(output.gpst,output.err_LTS(4,:),sz,output.err_LTS(4,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

%% Horizontal Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 1.5];
Marker = 'b.';

fig2 = figure(2); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.hor_err_LTS(1,:),sz,output.hor_err_LTS(1,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.hor_err_LTS(2,:),sz,output.hor_err_LTS(2,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(223); hold on; grid on
scatter(output.gpst,output.hor_err_LTS(3,:),sz,output.hor_err_LTS(3,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.hor_err_LTS(4,:),sz,output.hor_err_LTS(4,:), Marker)
ylim(errorLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')

%% GDOP Scatter Plot
sz = 15; % size of the dots
GDOPLimits = [0 inf];
Marker = 'b.';

fig3 = figure(3); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.GDOPLTS(1,:),sz,output.GDOPLTS(1,:), Marker)
ylim(GDOPLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(1))];
title(titleText,'Interpreter','tex')
ylabel('GDOP, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.GDOPLTS(2,:),sz,output.GDOPLTS(2,:), Marker)
ylim(GDOPLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(2))];
title(titleText,'Interpreter','tex')

subplot(223); hold on; grid on
scatter(output.gpst,output.GDOPLTS(3,:),sz,output.GDOPLTS(3,:), Marker)
ylim(GDOPLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(3))];
title(titleText,'Interpreter','tex')
ylabel('GDOP, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.GDOPLTS(4,:),sz,output.GDOPLTS(4,:), Marker)
ylim(GDOPLimits)
titleText = ['LTS Option = ',num2str(output.LTSOption(4))];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')
sgtitle('Geometric Dilution of Precision')



%% 3D Error CDF Plot
N = length(output.GDOPLS);
fig7 = figure(7); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','p','^','h','>','o'];
legend_list=cell(output.LTSn,1);
for idx = 1:output.LTSn
    [f,x] = ecdf(output.err_LTS(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('LTS Option: ', num2str(output.LTSOption(idx)));
end

leg7 = legend(legend_list); 
set(leg7,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of 3D Positioning Error');
xlabel('3D Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]); grid on
