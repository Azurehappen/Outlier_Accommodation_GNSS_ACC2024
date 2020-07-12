function log = save_result_linear(p,cpt,log,i,re_pos,clk_b,res)

log.pos_ecef{i} = re_pos;
ind_mark = cpt.svprn_mark ~= 0;
%------------------------%
[pos_llh,~,~]=ecef2llh_iter(re_pos.LS);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.Grdpos - re_pos.LS;
ned_err=R_e2g*err_pos;
log.ned_err_LS(i) = norm(ned_err);
log.hor_err_LS(i) = norm(ned_err(1:2));
log.err_LS(i) = norm(p.Grdpos - re_pos.LS);
log.GDOPLS(i) = cpt.GDOP.LS;
log.resLS(ind_mark,i) = res.LS;
%------------------------%
if p.eb_LSS == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.LSS);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.Grdpos - re_pos.LSS;
ned_err=R_e2g*err_pos;
log.ned_err_LSS(i) = norm(ned_err);
log.hor_err_LSS(i) = norm(ned_err(1:2));
log.err_LSS(i) = norm(p.Grdpos - re_pos.LSS);
log.GDOPLSS(i) = cpt.GDOP.LSS;
log.resLSS(ind_mark,i) = res.LSS;
end
%------------------------%
if p.eb_MShb == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.MShb);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.Grdpos - re_pos.MShb;
ned_err=R_e2g*err_pos;
log.ned_err_MShb(i) = norm(ned_err);
log.hor_err_MShb(i) = norm(ned_err(1:2));
log.err_MShb(i) = norm(p.Grdpos - re_pos.MShb);
log.GDOPMShb(i) = cpt.GDOP.MShb;
log.resMShb(ind_mark,i) = res.MShb;
end
%------------------------%
if p.eb_MStk == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.MStk);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.Grdpos - re_pos.MStk;
ned_err=R_e2g*err_pos;
log.ned_err_MStk(i) = norm(ned_err);
log.hor_err_MStk(i) = norm(ned_err(1:2));
log.err_MStk(i) = norm(p.Grdpos - re_pos.MStk);
log.GDOPMStk(i) = cpt.GDOP.MStk;
log.resMStk(ind_mark,i) = res.MStk;
end
%------------------------%
if p.eb_LTS == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.LTS);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.P_base - re_pos.LTS;
ned_err=R_e2g*err_pos;
log.ned_err_LTS(i) = norm(ned_err);
log.hor_err_LTS(i) = norm(ned_err(1:2));
log.err_LTS(i) = norm(p.P_base - re_pos.LTS);
log.GDOPLTS(i) = cpt.GDOP.LTS;
log.resLTS(ind_mark,i) = res.LTS;
end
%------------------------%
if p.eb_TD == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.TD1);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.P_base - re_pos.TD1;
ned_err=R_e2g*err_pos;
log.ned_err_TD1(i) = norm(ned_err);
log.hor_err_TD1(i) = norm(ned_err(1:2));
log.err_TD1(i) = norm(p.P_base - re_pos.TD1);
log.GDOPTD1(i) = cpt.GDOP.TD1;
log.resTD1(ind_mark,i) = res.TD1;
end
%------------------------%
if p.eb_TD == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.TD2);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.P_base - re_pos.TD2;
ned_err=R_e2g*err_pos;
log.ned_err_TD2(i) = norm(ned_err);
log.hor_err_TD2(i) = norm(ned_err(1:2));
log.err_TD2(i) = norm(p.P_base - re_pos.TD2);
log.GDOPTD2(i) = cpt.GDOP.TD2;
log.resTD2(ind_mark,i) = res.TD2;
end
%------------------------%
if p.eb_TD == 1
[pos_llh,~,~]=ecef2llh_iter(re_pos.TD3);
R_e2g=ll2R(pos_llh); % rotation matrix from ecef 2 geodetic frame
err_pos = p.P_base - re_pos.TD3;
ned_err=R_e2g*err_pos;
log.ned_err_TD3(i) = norm(ned_err);
log.hor_err_TD3(i) = norm(ned_err(1:2));
log.err_TD3(i) = norm(p.P_base - re_pos.TD3);
log.GDOPTD3(i) = cpt.GDOP.TD3;
log.resTD3(ind_mark,i) = res.TD3;
end
%------------------------%
log.rover_clk{i} = clk_b;
log.sv_num_GPS(i) = cpt.num_sv(1);log.sv_num_GLO(i) = cpt.num_sv(2);
log.sv_num_GAL(i) = cpt.num_sv(3);log.sv_num_BDS(i) = cpt.num_sv(4);

% start = 1; endi = log.num_obs_gps;
% log.res_GPS(:,i) = log.res(start:endi,i);
% start = start + log.num_obs_gps; endi = endi + log.num_obs_glo;
% log.res_GLO(:,i) = log.res(start:endi,i);
% start = start + log.num_obs_glo; endi = endi + log.num_obs_gal;
% log.res_GAL(:,i) = log.res(start:endi,i);
% start = start + log.num_obs_gal; endi = endi + log.num_obs_bds;
% log.res_BDS(:,i) = log.res(start:endi,i);

end