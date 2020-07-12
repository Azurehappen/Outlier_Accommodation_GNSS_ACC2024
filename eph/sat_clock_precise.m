function dt_sv_p = sat_clock_precise(p,IGSdata,prn,igs_idx,t_tsm)
% Compute precise satellite clock bias
% Input:
%       p: parameters
%       prn: SVID
%       igs_idx: the index in IGS data
%       t_tsm: time of transmit
dt_oc = limit_tgps(t_tsm - p.IGS.clk_iTOW(igs_idx)); 
if prn > size(IGSdata.clk_corr,1)
    dt_sv_p = 0;
else
    dt_sv_p=IGSdata.clk_corr(prn,igs_idx)+IGSdata.clk_vel(prn,igs_idx)*dt_oc+IGSdata.clk_acc(prn,igs_idx)*dt_oc^2; %%%% (eq.5) in meters
end
    dt_sv_p=dt_sv_p/p.c; %%% convert from meter to time 
end