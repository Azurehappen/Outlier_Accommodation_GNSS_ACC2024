function [cpt,n] = diff_corr_compute(p,cpt,obs_tr)
% Compute differential correction
% Input:
%       p       -- parameters which include base station's eph and obs
%       cpt     -- Satellite computation results
%       obs_tr  -- The time that rover received the data
diff_corr = NaN(length(cpt.corr_range),1);
ind_sv = find(cpt.svprn_mark~=0);
n = find(p.obs_b.tr_sow>(obs_tr-180) & p.obs_b.tr_sow<obs_tr);
if ~isempty(n)
%%%%%%%%%%%%%%%%%%%%%%%
delta = abs(p.obs_b.tr_sow(n)-obs_tr);
tridx = find(delta == min(delta));
tridx = n(tridx(end));
for i = 1:length(cpt.corr_range)
    prn = cpt.prn_record(ind_sv(i));
    sys_type = cpt.svprn_mark(ind_sv(i));
    obstr_B = p.obs_b.tr_sow(tridx);
    switch sys_type
        case 1
            pseR = p.obs_b.GPS(1).data.P(prn,tridx);
            Strength = p.obs_b.GPS(1).data.S(prn,tridx);
            message_duration = p.gps.message_duration;
            eph_info = p.eph_b.GPS;
        case 2
            pseR = p.obs_b.GLO(1).data.P(prn,tridx);
            Strength = p.obs_b.GLO(1).data.S(prn,tridx);
            message_duration = p.glo.message_duration;
            eph_info = p.eph_b.GLO;
            obstr_B = time_shift(obstr_B - p.glo.lps_gps);% Correct time diff from GPS time to GLO time
        case 3
            pseR = p.obs_b.GAL(1).data.P(prn,tridx);
            Strength = p.obs_b.GAL(1).data.S(prn,tridx);
            message_duration = p.gal.message_duration;
            eph_info = p.eph_b.GAL;
            obstr_B = time_shift(obstr_B - p.gal.lps_gps);% Correct time diff from GPS time to GAL time
        case 4
            pseR = p.obs_b.BDS(1).data.P(prn,tridx);
            Strength = p.obs_b.BDS(1).data.S(prn,tridx);
            message_duration = p.bds.message_duration;
            eph_info = p.eph_b.BDS;
            obstr_B = time_shift(obstr_B - p.bds.lps_gps);% Correct time diff from GPS time to BDS time
    end
    if (pseR~=0)&&(~isnan(pseR)) && prn<=size(eph_info.t_oc,1)
        tp_prime = pseR/p.c;
        t_sv = obstr_B - tp_prime;
        %---------------------------------% Find the time index in eph data
        teph = ephtidx(eph_info.t_oc{prn},t_sv,eph_info.SV_health(prn,:),message_duration);
        %---------------------------------%
        % Check the signal strength and sv health (health message o means ok)
        if ~isempty(teph) && Strength>=p.sig_strg
            switch(sys_type)
                % compute ephemeris to satellite position and clock bias
                case 1
                    [sat, dt_sv] = eph2pos(p,eph_info,p.obs_b,prn,teph,t_sv,'GPS');
                case 2
                    [sat, dt_sv] = geph2pos(p,eph_info,prn,teph,t_sv,'GLO');
                case 3
                    [sat, dt_sv] = eph2pos(p,eph_info,p.obs_b,prn,teph,t_sv,'GAL');
                case 4
                    [sat, dt_sv] = eph2pos(p,eph_info,p.obs_b,prn,teph,t_sv,'BDS'); 
            end
            range = norm(p.P_base-sat.pos_ecef)+sagnac(p,sat.pos_ecef,p.P_base);
            diff_corr(i) = pseR + p.c*dt_sv - range; % The correctoin here include receiver clock bias.
        else
            cpt.num_sv(sys_type) = cpt.num_sv(sys_type)-1;
            cpt.svprn_mark(ind_sv(i)) = 0;
            cpt.prn_record(ind_sv(i)) = 0;
        end
    else
        cpt.num_sv(sys_type) = cpt.num_sv(sys_type)-1;
        cpt.svprn_mark(ind_sv(i)) = 0;
        cpt.prn_record(ind_sv(i)) = 0;
    end
end
%%%%%%%%%%%%%%%%%%%%%%%
end
% Delete the data that has no diff correction
del_ind = find(isnan(diff_corr));
diff_corr(del_ind) = [];
cpt.corr_range(del_ind) = [];
cpt.s_pos_ecef(:,del_ind) = [];
cpt.s_v_ecef(:,del_ind) = [];
cpt.tp(del_ind) = [];
cpt.elev(del_ind) = [];
cpt.az(del_ind) = [];
cpt.sat_pos_Rcorr(:,del_ind) = [];
cpt.sat_posprc_Rcorr(:,del_ind) = [];
cpt.sat_v_Rcorr(:,del_ind) = [];
cpt.diff_corr = diff_corr;
end