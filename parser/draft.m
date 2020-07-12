log = log_temp.obs;
obs.tr_prime = [log.year';log.month';log.day';log.hour';log.min';log.sec'];
obs.tr_sow = log.iTOW';
obs.GAL=[];obs.GLO=[];obs.BDS=[];
obs.GPS.num_obstype = 4;
obs.GPS.f1 = 'L1';

len = size(log.sv_L1,1);
% P1 = zeros(32,len);
% 
% for i = 1:len
%     wh = log.sv_L1(i,:)~=0;
%     ind = log.sv_L1(i,wh);
%     P1(ind,i)=log.pr_L1(i,wh)';
% end
obs.GPS.S1 = log.snr_L1';
obs.GPS.P1 = log.pr_L1';
