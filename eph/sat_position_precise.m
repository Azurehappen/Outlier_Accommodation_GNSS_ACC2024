function sat_prec = sat_position_precise(p,IGSdata,sat_pos_ecef,sat_v_ecef,prn,igs_idx,t_tsm)
% Compute precise satellit postion
% Input:
%       p: parameters
%       prn: SVID
%       igs_idx: the index in IGS data
%       t_tsm: time of transmit


% orbit correction in radial along cross track direction
dP_rac = [IGSdata.orbit_x(prn,igs_idx);IGSdata.orbit_y(prn,igs_idx);IGSdata.orbit_z(prn,igs_idx)];
% orbit velocity correction in radial along cross track direction
dV_rac = [IGSdata.orbit_xv(prn,igs_idx);IGSdata.orbit_yv(prn,igs_idx);IGSdata.orbit_zv(prn,igs_idx)];

% theta = p.omge * limit_tgps(t_tsm - p.IGS.orbit_iTOW(igs_idx));
% R = [ cos(theta)  sin(theta)  0;
%      -sin(theta)  cos(theta)  0;
%       0               0       1];
% dP_rac = R*dP_rac;
% dV_rac = R*dV_rac;
% Compute position error
deph = limit_tgps(t_tsm - p.IGS.orbit_iTOW(igs_idx));
dP_rac = dP_rac + dV_rac*deph;

dP_ecef=RAC2ECEF(dP_rac,sat_pos_ecef,sat_v_ecef);

sat_prec = sat_pos_ecef-dP_ecef;

end