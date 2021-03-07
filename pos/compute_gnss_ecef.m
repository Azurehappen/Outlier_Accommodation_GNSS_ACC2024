function  log = compute_gnss_ecef(p,eph,obs)
% This function is to implement GNSS positioning with
% standard mode (without Iono, Trop, Es correction)
% or PPP mode.
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow); % The number of positioning points
% Initialize output
log.gpst = obs.tr_sow-obs.tr_sow(1);
log.err = NaN(1,N); % The position (Norm) error between estimated pos and true pos
log.hor_err = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
log.ned_err_norm = NaN(1,N); % NED frame norm error
log.ned_err = NaN(3,N); % NED frame error
log.pos_ecef = NaN(3,N); % Estimated position in ECEF
log.rover_clk = NaN(1,N); % Receiver clock bias (meter)
log.sv_num_GPS = NaN(1,N); % The amount of GPS satellite be used
log.sv_num_GLO = NaN(1,N); % The amount of GLO satellit  used
log.sv_num_GAL = NaN(1,N); % The amount of GAL satellite be used
log.sv_num_BDS = NaN(1,N); % The amount of BDS satellite be used
if ~isempty(obs.GPS)
    log.num_obs_gps = size(obs.GPS(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gps = 0;
end
if ~isempty(obs.GLO)
    log.num_obs_glo = size(obs.GLO(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_glo = 0;
end
if ~isempty(obs.GAL)
    log.num_obs_gal = size(obs.GAL(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gal = 0;
end
if ~isempty(obs.BDS)
    log.num_obs_bds = size(obs.BDS(1).data.P,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_bds = 0;
end
log.res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.elev_GPS = NaN(log.num_obs_gps,N); % The elevation of satellites
log.elev_GLO = NaN(log.num_obs_glo,N); 
log.elev_GAL = NaN(log.num_obs_gal,N); 
log.elev_BDS = NaN(log.num_obs_bds,N); 
log.res = [log.res_GPS;log.res_GAL;log.res_GLO;log.res_BDS];
log.elev = [log.elev_GPS;log.elev_GAL;log.elev_GLO;log.elev_BDS];
% Mark the sat prn that be computed
gpslog.svprn_mark = zeros(log.num_obs_gps,1);glolog.svprn_mark = zeros(log.num_obs_glo,1);
gallog.svprn_mark = zeros(log.num_obs_gal,1);bdslog.svprn_mark = zeros(log.num_obs_bds,1);
% Record the prn of each system for whose satllite been used
gpslog.prn_record = zeros(log.num_obs_gps,1);glolog.prn_record = zeros(log.num_obs_glo,1);
gallog.prn_record = zeros(log.num_obs_gal,1);bdslog.prn_record = zeros(log.num_obs_bds,1);
% Satellite position in ECEF frame
gpslog.s_pos_ecef = [];glolog.s_pos_ecef = [];gallog.s_pos_ecef = [];bdslog.s_pos_ecef = [];
if p.post_mode == 1
% precise satellite position in ECEF frame
gpslog.s_pos_prc = [];glolog.s_pos_prc = [];gallog.s_pos_prc = [];bdslog.s_pos_prc = [];
end
% Satellite velocity in ECEF frame
gpslog.s_v_ecef = [];glolog.s_v_ecef = [];gallog.s_v_ecef = [];bdslog.s_v_ecef = [];
% Signal propagation time
gpslog.tp = [];glolog.tp = [];gallog.tp = [];bdslog.tp = [];
% corrected pseudorange
gpslog.corr_range = [];glolog.corr_range = [];gallog.corr_range = [];bdslog.corr_range = [];
% The number of satellite be computed
gpslog.num_sv = 0;glolog.num_sv = 0;gallog.num_sv = 0;bdslog.num_sv = 0;
for i = 1:p.inval:N
    if ~isnan(p.Grdpos.t(1))
        index = round(p.Grdpos.t) == round(obs.tr_sow(i));
        grdpos = p.Grdpos.pos(:,index);
        if isempty(grdpos)
            continue;
        end
    else
        grdpos = p.Grdpos.pos;
    end
    if mod(i,4000)==0
        i
    end
    p.i = i; % To debug
    if ~isempty(obs.GPS)&& p.freq==1 && p.enableGPS == 1
        % GPS satellite position computation, Single frenquency receiver mode    
        gpslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gps,'GPS');
    end
    if ~isempty(obs.GLO)&& p.freq==1 && p.enableGLO == 1
        % GLO satellite position computation, Single frenquency receiver mode
        glolog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_glo,'GLO');
    end
    if ~isempty(obs.GAL)&& p.freq==1 && p.enableGAL == 1
        % GAL satellite position computation, Single frenquency receiver mode
        gallog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_gal,'GAL');
    end
    if ~isempty(obs.BDS)&& p.freq==1 && p.enableBDS == 1
        % BDS satellite position computation, Single frenquency receiver mode
        bdslog = satpost_corrpsedR_singlefreq(p,eph,obs,i,log.num_obs_bds,'BDS');
    end
    cpt.prn_record = [gpslog.prn_record;glolog.prn_record;gallog.prn_record;bdslog.prn_record];
    cpt.svprn_mark = [gpslog.svprn_mark;glolog.svprn_mark;gallog.svprn_mark;bdslog.svprn_mark];
    cpt.corr_range = [gpslog.corr_range;glolog.corr_range;gallog.corr_range;bdslog.corr_range];
    ind = find(cpt.corr_range==0);
    cpt.corr_range(ind) = [];
    cpt.s_pos_ecef = [gpslog.s_pos_ecef,glolog.s_pos_ecef,gallog.s_pos_ecef,bdslog.s_pos_ecef];
    cpt.s_pos_ecef(:,ind) = [];
    if p.post_mode == 1
        cpt.s_pos_prc = [gpslog.s_pos_prc,glolog.s_pos_prc,gallog.s_pos_prc,bdslog.s_pos_prc];
        cpt.s_pos_prc(:,ind) = [];
    end
    cpt.s_v_ecef = [gpslog.s_v_ecef,glolog.s_v_ecef,gallog.s_v_ecef,bdslog.s_v_ecef];
    cpt.s_v_ecef(:,ind) = [];
    cpt.tp = [gpslog.tp;glolog.tp;gallog.tp;bdslog.tp];
    cpt.tp(ind) = [];
    cpt.num_sv = [gpslog.num_sv,glolog.num_sv,gallog.num_sv,bdslog.num_sv];
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    if sum(cpt.num_sv)>=p.min_sv
        if p.state0(1)==0
            for kk = 1:3
                [re_pos,clock_bias,~] = userpos(p,cpt);
                p.state0 = [re_pos;clock_bias];
                if p.post_mode == 1
                    p.mk = 1;
                end
            end
        end
        % Rotate the sat pos to common reference frame
        cpt = earth_rotation_corr(p,cpt,p.state0(4)/p.c);
        % Check elevation
        cpt = elevaz_check(p,cpt,p.state0(1:3));
        if sum(cpt.num_sv)>=p.min_sv
            switch p.post_mode
                case 0 % Standard GNSS
%                     if p.elev_mark ==0
%                         log = save_result(p,cpt,log,i,re_pos,clock_bias,res);
%                     else
                        % Compute the final position
%                       [re_pos,clock_bias,res] = userpos_Rcorr(p,cpt);
                      tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                      [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                      cpt.IoFac = zeros(length(cpt.corr_range),1);
                      cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,[],rt);
                      cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                      [re_pos,clock_bias,res] = userpos(p,cpt);
                      [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
%                     end
                case 1 % PPP
                    tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                    [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                        %%-------------%%
                        cpt.IoFac = zeros(length(cpt.corr_range),1);
                        cpt = trop_iono_compute(p,eph,cpt,obs,p.state0(1:3),tdoy,p.USTEC,rt);
                        % Using the correction to the measurements
                        if ~isempty(find(cpt.iono_delay~=0, 1))
                        cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                        %%-------------%%
                        % Compute the final position
                        if p.double_diff == 0
                            [re_pos,clock_bias,res] = userpos(p,cpt);
                        else
                            [re_pos,clock_bias,res] = userpos_2diff(p,cpt);
                        end
                        if ~isempty(re_pos)
                            [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
                        end
                        end
                case 2 % DGNSS
                    if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                        [cpt,n] = diff_corr_compute(p,cpt,obs.tr_sow(i));
                        if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                            cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                            cpt.corr_range = cpt.corr_range - cpt.diff_corr;
%                             [re_pos,clock_bias,res] = userpos(p,cpt.sat_pos_Rcorr,...
%                                 cpt.corr_range,cpt.num_sv);
                            [re_pos,clock_bias,res] = userpos(p,cpt);
                            if ~isempty(re_pos)
                            [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
                            end
                        end
                    else
                        warning('No differential source given')
                    end
                    
                case 3 % VRS
                    tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                    [rt.week, rt.dow, rt.sow] = date2gnsst(obs.tr_prime(:,i)');
                    rt.sow = round(rt.sow);
                    cpt = vrs_corr_compute(p,cpt,eph,tdoy,rt);
                    if length(cpt.corr_range)>=p.min_sv
                        cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                        [re_pos,clock_bias,res] = userpos(p,cpt);
                        if ~isempty(re_pos)
                            [log,p.state0] = save_result(p,cpt,log,i,re_pos,clock_bias,res,grdpos);
                        end
                    end
                otherwise
                    warning('Unsupport positioning option');
            end
       end

    end
    
end






end