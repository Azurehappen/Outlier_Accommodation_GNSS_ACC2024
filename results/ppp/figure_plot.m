load('output_td_15.mat')
load('output_raps.mat')
load('output_ekf.mat')

% output_ekf = output;
% p_ekf = p;
% save('output_ekf.mat', 'output_ekf', 'p_ekf');
% output_td = output;
% p_td = p;
% save('output_td.mat', 'output_td', 'p_td');
% output_raps = output;
% p_raps = p;
% save('output_raps.mat', 'output_raps', 'p_raps');
%%
set(0,'defaultfigurecolor','w')
figure
plot(p_raps.t,output_raps.raps_flag)
title('RAPS flag: 0 means not used, TD instead')
xlabel('Time')
ylabel('Flag')
grid on

ind = output_raps.raps_flag == 1;
raps_hor_err = output_raps.hor_err(ind);
raps_ver_err = abs(output_raps.ned_err(3,ind));
raps_cost = output_raps.cost(ind);
td_hor_err = output_td.hor_err(ind);
td_ver_err = abs(output_td.ned_err(3,ind));
td_cost = output_td.cost(ind);
ekf_hor_err = output_ekf.hor_err(ind);
ekf_ver_err = abs(output_ekf.ned_err(3,ind));
ekf_cost = output_ekf.cost(ind);

total_ekf = output_ekf.num_meas_used;
total_td = output_td.num_meas_used;
total_raps = output_raps.num_meas_used;
% total_td = total_ekf - total_td;
% total_raps = total_ekf - total_raps;

total_ekf = total_ekf(ind);
total_td = total_td(ind);
total_raps = total_raps(ind);

time = p_raps.t;
st_t = time(330);
ed_t = time(1140);
time = time(ind);
st_t_i = find(time == st_t);
ed_t_i = find(time == ed_t);

disp('Overall performance')
nonNaNCount = sum(~isnan(raps_hor_err));
fprintf('RAPS Hor <= 1.0 m: %.2f%%\n', sum(raps_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Hor <= 1.5 m: %.2f%%\n', sum(raps_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Ver <= 3.0 m: %.2f%%\n', sum(raps_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('RAPS Hor Mean: %.2f\n', mean(raps_hor_err));
fprintf('RAPS Hor RMS: %.2f\n', rms(raps_hor_err));
fprintf('RAPS Hor Max: %.2f\n', max(raps_hor_err)); 
fprintf('RAPS Ver Mean: %.2f\n', mean(raps_ver_err));
fprintf('RAPS Ver RMS: %.2f\n', rms(raps_ver_err));
fprintf('RAPS Ver Max: %.2f\n', max(raps_ver_err));

nonNaNCount = sum(~isnan(td_hor_err));
fprintf('TD Hor <= 1.0 m: %.2f%%\n', sum(td_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('TD Hor <= 1.5 m: %.2f%%\n', sum(td_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('TD Ver <= 3.0 m: %.2f%%\n', sum(td_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('TD Hor Mean: %.2f\n', mean(td_hor_err));
fprintf('TD Hor RMS: %.2f\n', rms(td_hor_err));
fprintf('TD Hor Max: %.2f\n', max(td_hor_err));
fprintf('TD Ver Mean: %.2f\n', mean(td_ver_err));
fprintf('TD Ver RMS: %.2f\n', rms(td_ver_err));
fprintf('TD Ver Max: %.2f\n', max(td_ver_err));

nonNaNCount = sum(~isnan(ekf_hor_err));
fprintf('EKF Hor <= 1.0 m: %.2f%%\n', sum(ekf_hor_err <= 1.0) / nonNaNCount * 100);
fprintf('EKF Hor <= 1.5 m: %.2f%%\n', sum(ekf_hor_err <= 1.5) / nonNaNCount * 100);
fprintf('EKF Ver <= 3.0 m: %.2f%%\n', sum(ekf_ver_err <= 3.0) / nonNaNCount * 100);
fprintf('EKF Hor Mean: %.2f\n', mean(ekf_hor_err));
fprintf('EKF Hor RMS: %.2f\n', rms(ekf_hor_err));
fprintf('EKF Hor Max: %.2f\n', max(ekf_hor_err));
fprintf('EKF Ver Mean: %.2f\n', mean(ekf_ver_err));
fprintf('EKF Ver RMS: %.2f\n', rms(ekf_ver_err));
fprintf('EKF Ver Max: %.2f\n', max(ekf_ver_err));

% Split for open-sky or non-open-sky
% Non-open-sky: ns
i_ns = 333:1140;
raps_flag_ns = output_raps.raps_flag(i_ns);
raps_hor_err_ns = output_raps.hor_err(i_ns);
raps_ver_err_ns = abs(output_raps.ned_err(3,i_ns));
td_hor_err_ns = output_td.hor_err(i_ns);
td_ver_err_ns = abs(output_td.ned_err(3,i_ns));
ekf_hor_err_ns = output_ekf.hor_err(i_ns);
ekf_ver_err_ns = abs(output_ekf.ned_err(3,i_ns));

ind = raps_flag_ns == 1;
raps_hor_err_ns = raps_hor_err_ns(ind);
raps_ver_err_ns = raps_ver_err_ns(ind);
td_hor_err_ns = td_hor_err_ns(ind);
td_ver_err_ns = td_ver_err_ns(ind);
ekf_hor_err_ns = ekf_hor_err_ns(ind);
ekf_ver_err_ns = ekf_ver_err_ns(ind);

disp('Non-open sky Overall performance')
nonNaNCount = sum(~isnan(raps_hor_err_ns));
fprintf('RAPS Hor <= 1.0 m: %.2f%%\n', sum(raps_hor_err_ns <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Hor <= 1.5 m: %.2f%%\n', sum(raps_hor_err_ns <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Ver <= 3.0 m: %.2f%%\n', sum(raps_ver_err_ns <= 3.0) / nonNaNCount * 100);
fprintf('RAPS Hor Mean: %.2f\n', mean(raps_hor_err_ns));
fprintf('RAPS Hor RMS: %.2f\n', rms(raps_hor_err_ns));
fprintf('RAPS Hor Max: %.2f\n', max(raps_hor_err_ns));
fprintf('RAPS Ver Mean: %.2f\n', mean(raps_ver_err_ns));
fprintf('RAPS Ver RMS: %.2f\n', rms(raps_ver_err_ns));
fprintf('RAPS Ver Max: %.2f\n', max(raps_ver_err_ns));

nonNaNCount = sum(~isnan(td_hor_err_ns));
fprintf('TD Hor <= 1.0 m: %.2f%%\n', sum(td_hor_err_ns <= 1.0) / nonNaNCount * 100);
fprintf('TD Hor <= 1.5 m: %.2f%%\n', sum(td_hor_err_ns <= 1.5) / nonNaNCount * 100);
fprintf('TD Ver <= 3.0 m: %.2f%%\n', sum(td_ver_err_ns <= 3.0) / nonNaNCount * 100);
fprintf('TD Hor Mean: %.2f\n', mean(td_hor_err_ns));
fprintf('TD Hor RMS: %.2f\n', rms(td_hor_err_ns));
fprintf('TD Hor Max: %.2f\n', max(td_hor_err_ns));
fprintf('TD Ver Mean: %.2f\n', mean(td_ver_err_ns));
fprintf('TD Ver RMS: %.2f\n', rms(td_ver_err_ns));
fprintf('TD Ver Max: %.2f\n', max(td_ver_err_ns));

nonNaNCount = sum(~isnan(ekf_hor_err_ns));
fprintf('EKF Hor <= 1.0 m: %.2f%%\n', sum(ekf_hor_err_ns <= 1.0) / nonNaNCount * 100);
fprintf('EKF Hor <= 1.5 m: %.2f%%\n', sum(ekf_hor_err_ns <= 1.5) / nonNaNCount * 100);
fprintf('EKF Ver <= 3.0 m: %.2f%%\n', sum(ekf_ver_err_ns <= 3.0) / nonNaNCount * 100);
fprintf('EKF Hor Mean: %.2f\n', mean(ekf_hor_err_ns));
fprintf('EKF Hor RMS: %.2f\n', rms(ekf_hor_err_ns));
fprintf('EKF Hor Max: %.2f\n', max(ekf_hor_err_ns));
fprintf('EKF Ver Mean: %.2f\n', mean(ekf_ver_err_ns));
fprintf('EKF Ver RMS: %.2f\n', rms(ekf_ver_err_ns));
fprintf('EKF Ver Max: %.2f\n', max(ekf_ver_err_ns));

% Open-sky: os
i_os = [1:332,1141:length(output_raps.raps_flag)];
raps_flag_os = output_raps.raps_flag(i_os);
raps_hor_err_os = output_raps.hor_err(i_os);
raps_ver_err_os = abs(output_raps.ned_err(3,i_os));
td_hor_err_os = output_td.hor_err(i_os);
td_ver_err_os = abs(output_td.ned_err(3,i_os));
ekf_hor_err_os = output_ekf.hor_err(i_os);
ekf_ver_err_os = abs(output_ekf.ned_err(3,i_os));

ind = raps_flag_os == 1;
raps_hor_err_os = raps_hor_err_os(ind);
raps_ver_err_os = raps_ver_err_os(ind);
td_hor_err_os = td_hor_err_os(ind);
td_ver_err_os = td_ver_err_os(ind);
ekf_hor_err_os = ekf_hor_err_os(ind);
ekf_ver_err_os = ekf_ver_err_os(ind);

disp('Open sky Overall performance')
nonNaNCount = sum(~isnan(raps_hor_err_os));
fprintf('RAPS Hor <= 1.0 m: %.2f%%\n', sum(raps_hor_err_os <= 1.0) / nonNaNCount * 100);
fprintf('RAPS Hor <= 1.5 m: %.2f%%\n', sum(raps_hor_err_os <= 1.5) / nonNaNCount * 100);
fprintf('RAPS Ver <= 3.0 m: %.2f%%\n', sum(raps_ver_err_os <= 3.0) / nonNaNCount * 100);
fprintf('RAPS Hor Mean: %.2f\n', mean(raps_hor_err_os));
fprintf('RAPS Hor RMS: %.2f\n', rms(raps_hor_err_os));
fprintf('RAPS Hor Max: %.2f\n', max(raps_hor_err_os));
fprintf('RAPS Ver Mean: %.2f\n', mean(raps_ver_err_os));
fprintf('RAPS Ver RMS: %.2f\n', rms(raps_ver_err_os));
fprintf('RAPS Ver Max: %.2f\n', max(raps_ver_err_os));

nonNaNCount = sum(~isnan(td_hor_err_os));
fprintf('TD Hor <= 1.0 m: %.2f%%\n', sum(td_hor_err_os <= 1.0) / nonNaNCount * 100);
fprintf('TD Hor <= 1.5 m: %.2f%%\n', sum(td_hor_err_os <= 1.5) / nonNaNCount * 100);
fprintf('TD Ver <= 3.0 m: %.2f%%\n', sum(td_ver_err_os <= 3.0) / nonNaNCount * 100);
fprintf('TD Hor Mean: %.2f\n', mean(td_hor_err_os));
fprintf('TD Hor RMS: %.2f\n', rms(td_hor_err_os));
fprintf('TD Hor Max: %.2f\n', max(td_hor_err_os));
fprintf('TD Ver Mean: %.2f\n', mean(td_ver_err_os));
fprintf('TD Ver RMS: %.2f\n', rms(td_ver_err_os));
fprintf('TD Ver Max: %.2f\n', max(td_ver_err_os));

nonNaNCount = sum(~isnan(ekf_hor_err_os));
fprintf('EKF Hor <= 1.0 m: %.2f%%\n', sum(ekf_hor_err_os <= 1.0) / nonNaNCount * 100);
fprintf('EKF Hor <= 1.5 m: %.2f%%\n', sum(ekf_hor_err_os <= 1.5) / nonNaNCount * 100);
fprintf('EKF Ver <= 3.0 m: %.2f%%\n', sum(ekf_ver_err_os <= 3.0) / nonNaNCount * 100);
fprintf('EKF Hor Mean: %.2f\n', mean(ekf_hor_err_os));
fprintf('EKF Hor RMS: %.2f\n', rms(ekf_hor_err_os));
fprintf('EKF Hor Max: %.2f\n', max(ekf_hor_err_os));
fprintf('EKF Ver Mean: %.2f\n', mean(ekf_ver_err_os));
fprintf('EKF Ver RMS: %.2f\n', rms(ekf_ver_err_os));
fprintf('EKF Ver Max: %.2f\n', max(ekf_ver_err_os));
%%
purple = [0.4940, 0.1840, 0.5560];
blue = [0, 0.4470, 0.7410];
red = [0.6350, 0.0780, 0.1840];
green = [0.4660, 0.6740, 0.1880];
orange = [0.9290, 0.6940, 0.1250];

ekf_color = blue;
td_color = orange;
raps_color = red;

epoch = 1:length(raps_cost);
% Plot cost
figure(1)
subplot(211)
plot(epoch, ekf_cost, '.', 'Color', ekf_color)
hold on
plot(epoch, td_cost, '.', 'Color', td_color)
plot(epoch, raps_cost, '.', 'Color', raps_color)
% Add a green box
yMin = 0.01;
yMax = 10^4;
rectangle('Position', [epoch(st_t_i), yMin, epoch(ed_t_i)-epoch(st_t_i), yMax-yMin], 'EdgeColor', 'g')
hold off
set(gca, 'XGrid', 'off', 'YGrid', 'on')
axis tight
legend('EKF', 'TD', 'RAPS')
ylabel('Risk')
xlabel('Epochs, sec')
set(gca, 'YScale', 'log')

subplot(212)
yyaxis left
plot(epoch, total_ekf-total_td, '.', 'Color', td_color)
hold on
plot(epoch, total_ekf-total_raps, '.', 'Color', raps_color)
ylabel('No. of Meas. being removed')
set(gca, 'YColor', 'k')
% Right y-axis for total measurements
yyaxis right
plot(epoch, total_ekf, 'Color', ekf_color)
set(gca, 'YColor', ekf_color)
ylabel('Total Meas.')
hold off
grid on
axis tight
legend('TD', 'RAPS', 'Total')
xlabel('Epochs, sec')

% Plot pos errors
epoch = 1:length(ekf_hor_err);
figure(2)
subplot(211)
plot(epoch, ekf_hor_err, "Color", ekf_color)
hold on
plot(epoch, td_hor_err, "Color", td_color)
plot(epoch, raps_hor_err, "Color", raps_color)
hold off
grid on
axis tight
legend('EKF', 'TD', 'RAPS')
ylabel('Horizontal errors, meters')
xlabel('Epochs, sec')

subplot(212)
plot(epoch, ekf_ver_err, "Color", ekf_color)
hold on
plot(epoch, td_ver_err, "Color", td_color)
plot(epoch, raps_ver_err, "Color", raps_color)
hold off
grid on
axis tight
legend('EKF', 'TD', 'RAPS')
ylabel('Vertical errors, meters')
xlabel('Epochs, sec')

L = length(ekf_hor_err);
[f_ekf,x_ekf] = ecdf(ekf_hor_err);
[f_td,x_td] = ecdf(td_hor_err);
[f_raps,x_raps] = ecdf(raps_hor_err);
figure(3)
subplot(211)
h_ekf = semilogx(x_ekf,f_ekf); hold on;
h_ekf.Color = ekf_color;
h_ekf.Marker = '*';
h_ekf.LineWidth = 0.2;
h_ekf.MarkerIndices = 1:200:L;

h_td = semilogx(x_td,f_td); hold on;
h_td.Color = td_color;
h_td.Marker = 'o';
h_td.LineWidth = 0.2;
h_td.MarkerIndices = 1:200:L;

h_raps = semilogx(x_raps,f_raps); hold on;
h_raps.Color = raps_color;
h_raps.Marker = '>';
h_raps.LineWidth = 0.2;
h_raps.MarkerIndices = 1:200:L;
xline(1.0)
xline(1.5)
grid on
xlim([0.01, 10])
legend([h_ekf,h_td,h_raps],{'EKF', 'TD','RAPS'})
legend('location','best');
currentXTicks = get(gca, 'XTick');
set(gca, 'XTick', unique([currentXTicks, 1.5]));
%title('CDF of Horizontal Error');
xlabel('Horizontal error, meter');
ylabel('Cumulative Probability');

[f_ekf,x_ekf] = ecdf(ekf_ver_err);
[f_td,x_td] = ecdf(td_ver_err);
[f_raps,x_raps] = ecdf(raps_ver_err);
subplot(212)
h_ekf = semilogx(x_ekf,f_ekf); hold on;
h_ekf.Color = ekf_color;
h_ekf.Marker = '*';
h_ekf.LineWidth = 0.2;
h_ekf.MarkerIndices = 1:200:L;

h_td = semilogx(x_td,f_td); hold on;
h_td.Color = td_color;
h_td.Marker = 'o';
h_td.LineWidth = 0.2;
h_td.MarkerIndices = 1:200:L;

h_raps = semilogx(x_raps,f_raps); hold on;
h_raps.Color = raps_color;
h_raps.Marker = '>';
h_raps.LineWidth = 0.2;
h_raps.MarkerIndices = 1:200:L;
xline(3)
grid on
legend([h_ekf,h_td,h_raps],{'EKF', 'TD','RAPS'})
legend('location','best');
xlim([0.01, 10])
currentXTicks = get(gca, 'XTick');
set(gca, 'XTick', unique([currentXTicks, 3.0]));
%title('CDF of Horizontal Error');
xlabel('Vertical error, meter');
ylabel('Cumulative Probability');
%%
% Postiro measurements residual
% Prior measuremnts residual

