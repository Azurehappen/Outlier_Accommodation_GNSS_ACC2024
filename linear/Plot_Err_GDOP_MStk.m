%% Figures
load('./../output_linear.mat')

%% 3D Error CDF Plot
N = length(output.GDOPLS);
fig7 = figure(8); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','p','^','h','>','o'];
legend_list=cell(output.MStkn,1);
for idx = 1:output.MStkn
    [f,x] = ecdf(output.err_MStk(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('MStk BisqConst = ', num2str(output.MStkConst(idx)));
end

leg7 = legend(legend_list); 
set(leg7,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of 3D Positioning Error');
xlabel('3D Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]); grid on

%% 3D Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 4];
Marker = 'b.';

fig1 = figure(1); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.err_MStk(1,:),sz,output.err_MStk(1,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator (Tukey), BisqConst = ',num2str(output.MStkConst(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.err_MStk(2,:),sz,output.err_MStk(2,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('3D positioning error')

subplot(223); hold on; grid on
scatter(output.gpst,output.err_MStk(3,:),sz,output.err_MStk(3,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator (Tukey), BisqConst = ',num2str(output.MStkConst(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(224); hold on; grid on
scatter(output.gpst,output.err_MStk(4,:),sz,output.err_MStk(4,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator (Tukey), BisqConst = ',num2str(output.MStkConst(4))];
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
scatter(output.gpst,output.GDOPMStk(1,:),sz,output.GDOPMStk(1,:), Marker)
ylim(GDOPLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(1))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.GDOPMStk(2,:),sz,output.GDOPMStk(2,:), Marker)
ylim(GDOPLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(2))];
title(titleText,'Interpreter','tex')

subplot(223); hold on; grid on
scatter(output.gpst,output.GDOPMStk(3,:),sz,output.GDOPMStk(3,:), Marker)
ylim(GDOPLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(3))];
title(titleText,'Interpreter','tex')
ylabel('DOP, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.GDOPMStk(4,:),sz,output.GDOPMStk(4,:), Marker)
ylim(GDOPLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(4))];
title(titleText,'Interpreter','tex')
xlabel('Receiver time using GPS second')
sgtitle('GDOP')

%% Horizontal Error CDF Plot
N = length(output.GDOPLS);
fig6 = figure(6); clf
colorlist = ['r','g','b','c','m','k','r','g','b','c','m','k'];
markerlist = ['o','>','d','v','s','<','x','>','p','h','^','o'];
legend_list=cell(output.MStkn,1);
for idx = 1:output.MStkn
    [f,x] = ecdf(output.hor_err_MStk(idx,:));
    h_mer = semilogx(x,f); hold on;
    h_mer.Color = colorlist(idx);
    h_mer.Marker = markerlist(idx);
    h_mer.LineWidth = 1;
    h_mer.MarkerIndices = 1:5e3:N; %1:200:N;
    legend_list{idx}=strcat('MStk BisqConst = ', num2str(output.MStkConst(idx)));
end

leg6 = legend(legend_list); 
set(leg6,'Location','northeastoutside');
title('')
sgtitle('Cumulative Distribution Of Horizontal Positioning Error');
xlabel('Horizontal Positioning Error, m');
ylabel('Cumulative Probability');
xlim([1e-2 10]); grid on

%% Horizontal Positioning error Scatter Plot
sz = 15; % size of the dots
errorLimits = [0 1.5];
Marker = 'b.';

fig2 = figure(2); clf
colormap jet
subplot(221); hold on; grid on
scatter(output.gpst,output.hor_err_MStk(1,:),sz,output.hor_err_MStk(1,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(1))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(222); hold on; grid on
scatter(output.gpst,output.hor_err_MStk(2,:),sz,output.hor_err_MStk(2,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(2))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')

subplot(223); hold on; grid on
scatter(output.gpst,output.hor_err_MStk(3,:),sz,output.hor_err_MStk(3,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(3))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')

subplot(224); hold on; grid on
scatter(output.gpst,output.hor_err_MStk(4,:),sz,output.hor_err_MStk(4,:), Marker)
ylim(errorLimits)
titleText = ['M-Estimator(Tukey), BisqConst = ',num2str(output.MStkConst(4))];
title(titleText,'Interpreter','tex')
ylabel('Error, unit: meters')
xlabel('Receiver time using GPS second')
sgtitle('Horizontal positioning error')


%%
err_MStk = zeros(output.MStkn,1);
for idx = 1:output.MStkn
    err_MStk = mean(output.err_MStk(idx,:));    
end

err_MShb = zeros(output.MShbn,1);
for idx = 1:output.MShbn
    err_MShb = mean(output.err_MShb(idx,:));    
end

err_MStk = zeros(output.MStkn,1);
for idx = 1:output.MStkn
    err_MStk = mean(output.err_MStk(idx,:));    
end

err_MStk = zeros(output.MStkn,1);
for idx = 1:output.MStkn
    err_MStk = mean(output.err_MStk(idx,:));    
end

err_MStk = zeros(output.MStkn,1);
for idx = 1:output.MStkn
    err_MStk = mean(output.err_MStk(idx,:));    
end
