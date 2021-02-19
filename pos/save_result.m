function [log,state] = save_result(p,cpt,log,i,re_pos,clk_b,res,grdpos)

log.pos_ecef(:,i) = re_pos;
state = [re_pos;clk_b];
%------------------------%
[pos_llh,~,~]=ecef2llh_iter(re_pos);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = grdpos - re_pos;
ned_err=R_e2g*err_pos;
log.ned_err(:,i) = ned_err;
%------------------------%
log.ned_err_norm(i) = norm(ned_err);
log.hor_err(i) = norm(ned_err(1:2));
log.err(i) = norm(grdpos - re_pos);
log.rover_clk(i) = clk_b;
log.sv_num_GPS(i) = cpt.num_sv(1);log.sv_num_GLO(i) = cpt.num_sv(2);
log.sv_num_GAL(i) = cpt.num_sv(3);log.sv_num_BDS(i) = cpt.num_sv(4);
ind_mark = cpt.svprn_mark ~= 0;
log.res(ind_mark,i) = res;
log.elev(ind_mark,i) = cpt.elev;
start = 1; endi = log.num_obs_gps;
log.res_GPS(:,i) = log.res(start:endi,i);
log.elev_GPS(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gps; endi = endi + log.num_obs_glo;
log.res_GLO(:,i) = log.res(start:endi,i);
log.elev_GLO(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_glo; endi = endi + log.num_obs_gal;
log.res_GAL(:,i) = log.res(start:endi,i);
log.elev_GAL(:,i) = log.elev(start:endi,i);
start = start + log.num_obs_gal; endi = endi + log.num_obs_bds;
log.res_BDS(:,i) = log.res(start:endi,i);
log.elev_BDS(:,i) = log.elev(start:endi,i);

end