% SCRIPT FOR PLOTTING OUTPUTS OF DIFFERENT METHODS 
% This is s that we don't have to do computations every time we need
%   to plot things. -Bernard 
% ===================================================
clear; close all; 

% For Debug Purposes
MER = 1; 

% load outputs
LSoutput = load('./linear/output/LSoutput.mat'); 
LTSoutput = load('./linear/output/LTSoutput.mat'); 
MERoutput= load('./linear/output/MERoutput.mat'); 

%% PLOT
figure(1)
output = LSoutput.LSoutput; 
scatter(output.gpst,output.hor_err,'b.')
title('LS Horizontal Positioning Error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');
grid on; hold on

figure(2)
scatter(output.gpst,output.sv_num_GPS+output.sv_num_GLO+output.sv_num_GAL+output.sv_num_BDS,'b.')
title('Total Satellites Used')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');
grid on; hold on

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

figure(3)
plot(output.gpst,output.cost,'b.')
title('Measurement Residual Comparison')
xlabel('Receiver time using GPS second')
title('Cost vs time')
xlabel('Receiver time using GPS second')
ylabel('Cost: $$\sum_{i=1}^{(m+n)} (\textbf{r}_i(\textbf{x}))^2$$',...
    'Interpreter','Latex')
grid on; hold on


%%
output = LTSoutput.LTSoutput; 
figure(1)
scatter(output.gpst,output.hor_err,'r.')
title('LTS Horizontal positioning error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');
grid on; hold on

figure(2)
scatter(output.gpst,output.sv_num_GPS+output.sv_num_GLO+output.sv_num_GAL+output.sv_num_BDS,'r.')
title('Total Satellites Used')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');
grid on; hold on

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

figure(3)
plot(output.gpst,output.cost,'r.')
title('Measurement residual comparison')
xlabel('Receiver time using GPS second')
title('Cost vs time')
xlabel('Receiver time using GPS second')
ylabel('Cost: $$\sum_{i=1}^{(m+n)} (\textbf{r}_i(\textbf{x}))^2$$',...
    'Interpreter','Latex'); 
grid on; hold on 


%% 
MER =1; 
if MER 
    
output = MERoutput.output; 
figure(1)
scatter(output.gpst,output.hor_err,'g.')
title('Horizontal positioning error')
xlabel('Receiver time using GPS second')
ylabel('Error, unit: meter');
grid on; hold on

figure(2)
scatter(output.gpst,output.sv_num_GPS+output.sv_num_GLO+output.sv_num_GAL+output.sv_num_BDS,'g.')
title('Total Satellites Used')
xlabel('Receiver time using GPS second')
ylabel('Distance, unit: meter');
grid on; hold on

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

figure(3)
plot(output.gpst,output.cost,'g.')
title('Measurement residual comparison')
xlabel('Receiver time using GPS second')
title('Cost vs time')
xlabel('Receiver time using GPS second')
ylabel('Cost: $$\sum_{i=1}^{(m+n)} (\textbf{r}_i(\textbf{x}))^2$$',...
    'Interpreter','Latex'); 
grid on; hold on 

end 
    
%%
for i = 1:length(cost)-1
cost_diff(i)= cost(i+1)-cost(i); 

end 