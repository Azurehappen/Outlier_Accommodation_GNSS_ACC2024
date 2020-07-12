function  log = linear_gnss_ecef(p,eph,obs)
% This function is to implement GNSS positioning with
% linear mode (without Iono, Trop, Es correction)
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow); % The number of positioning points
% Initialize output
log.gpst = obs.tr_sow-obs.tr_sow(1);
log.err_LS = NaN(1,N); % The position (Norm) error between estimated pos and true pos
log.hor_err_LS = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
log.ned_err_LS = NaN(1,N); % NED frame error
log.err_LSS = NaN(1,N); log.hor_err_LSS = NaN(1,N); log.ned_err_LSS = NaN(1,N);
log.err_MShb = NaN(1,N); log.hor_err_MShb = NaN(1,N); log.ned_err_MShb = NaN(1,N);
log.err_MStk = NaN(1,N); log.hor_err_MStk = NaN(1,N); log.ned_err_MStk = NaN(1,N);
log.err_LTS = NaN(1,N); log.hor_err_LTS = NaN(1,N); log.ned_err_LTS = NaN(1,N);
log.err_TD1 = NaN(1,N); log.hor_err_TD1 = NaN(1,N); log.ned_err_TD1 = NaN(1,N);
log.err_TD2 = NaN(1,N); log.hor_err_TD2 = NaN(1,N); log.ned_err_TD2 = NaN(1,N);
log.err_TD3 = NaN(1,N); log.hor_err_TD3 = NaN(1,N); log.ned_err_TD3 = NaN(1,N);
log.pos_ecef = cell(1,N); % Estimated position in ECEF
log.rover_clk = cell(1,N); % Receiver clock bias (meter)
log.GDOPLS = NaN(1,N); % GDOP function (meter)
log.GDOPLSS = NaN(1,N);
log.GDOPMShb = NaN(1,N);
log.GDOPMStk = NaN(1,N);
log.GDOPLTS = NaN(1,N);
log.GDOPTD1 = NaN(1,N);
log.GDOPTD2 = NaN(1,N);
log.GDOPTD3 = NaN(1,N);
log.sv_num_GPS = NaN(1,N); % The amount of GPS satellite be used
log.sv_num_GLO = NaN(1,N); % The amount of GLO satellite be used
log.sv_num_GAL = NaN(1,N); % The amount of GAL satellite be used
log.sv_num_BDS = NaN(1,N); % The amount of BDS satellite be used

if ~isempty(obs.GPS)
    log.num_obs_gps = size(obs.GPS.P1,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gps = 0;
end
if ~isempty(obs.GLO)
    log.num_obs_glo = size(obs.GLO.P1,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_glo = 0;
end
if ~isempty(obs.GAL)
    log.num_obs_gal = size(obs.GAL.P1,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_gal = 0;
end
if ~isempty(obs.BDS)
    log.num_obs_bds = size(obs.BDS.P1,1); % The maximum of PRN recorded in obs data
else
    log.num_obs_bds = 0;
end
log.res_GPS = NaN(log.num_obs_gps,N); % The residual at the end
log.res_GLO = NaN(log.num_obs_glo,N); % The residual at the end
log.res_GAL = NaN(log.num_obs_gal,N); % The residual at the end
log.res_BDS = NaN(log.num_obs_bds,N); % The residual at the end
log.resLS = [log.res_GPS;
log.res_GAL;
log.res_GLO;log.res_BDS];
log.resLSS = log.resLS;
log.resMShb = log.resLS;
log.resMStk = log.resLS;
log.resLTS = log.resLS;
log.resTD1 = log.resLS;
log.resTD2 = log.resLS;
log.resTD3 = log.resLS;
% Mark the sat prn that be computed
gpslog.svprn_mark = [];
glolog.svprn_mark = [];
gallog.svprn_mark = [];
bdslog.svprn_mark = [];
% Record the prn of each system for whose satllite been used
gpslog.prn_record = [];
glolog.prn_record = [];
gallog.prn_record = [];
bdslog.prn_record = [];
% Satellite position in ECEF frame
gpslog.s_pos_ecef = [];
glolog.s_pos_ecef = [];
gallog.s_pos_ecef = [];
bdslog.s_pos_ecef = [];
% corrected pseudorange
gpslog.corr_range = [];
glolog.corr_range = [];
gallog.corr_range = [];
bdslog.corr_range = [];
% The number of satellite be computed
gpslog.num_sv = 0;
glolog.num_sv = 0;
gallog.num_sv = 0;
bdslog.num_sv = 0;
for i = 1:1000
    p.i = i;
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
    cpt.num_sv = [gpslog.num_sv,glolog.num_sv,gallog.num_sv,bdslog.num_sv];
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    
    if sum(cpt.num_sv)>=p.min_sv
%         if p.state0(1)==0
%             [A,p.x0(4),~] = userpos(p,cpt.s_pos_ecef,cpt.corr_range,cpt.num_sv);
%         end
        % Check elevation
        re_pos = p.x0(1:3);
        cpt = elevaz_check_linear(p,cpt,re_pos);
        if sum(cpt.num_sv)>=p.min_sv
            switch p.post_mode
                case 0 % Standard GNSS
                    if p.elev_mark ==0
                        log = save_result_linear(p,cpt,log,i,re_pos,clock_bias,res);
                    else
                        % Compute the final position
                        [re_pos,clock_bias,res,cpt.GDOP] = userposlinear(p,cpt.s_pos_ecef,...
                            cpt.corr_range,cpt.num_sv);
                        log = save_result_linear(p,cpt,log,i,re_pos,clock_bias,res);
                    end
                case 1 % PPP
                    tdoy = doy(obs.tr_prime(1:3,i)); % Day of year
                    user_time = obs.tr_prime(3:6,i)';
                    %%-------------%%
                    cpt = trop_iono_compute(p,cpt,obs,re_pos,tdoy,p.USTEC,user_time);
                    % Using the correction to the measurements
                    cpt.corr_range = cpt.corr_range - cpt.trop_delay - cpt.iono_delay;
                    %%-------------%%
                    % Compute the final position
                    [re_pos,clock_bias,res,cpt.GDOP] = userposlinear(p,cpt.s_pos_ecef,...
                        cpt.corr_range,cpt.num_sv);
                    log = save_result_linear(p,cpt,log,i,re_pos,clock_bias,res);
                case 2 % DGNSS
                    if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                        [cpt,n] = diff_corr_compute(p,cpt,obs.tr_sow(i));
                        if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                            cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                            cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                            %---------% Get prior for reciver clock bias
                            [~,p.x0(4),~] = userpos(p,cpt);
                            %---------%
                            [re_pos,clock_bias,res,cpt.GDOP] = userposlinear(p,cpt.s_pos_ecef,...
                                cpt.corr_range,cpt.num_sv);
                            log = save_result_linear(p,cpt,log,i,re_pos,clock_bias,res);
                        end
                    else
                        warning('No differential source given')
                    end
                otherwise
                    warning('Unsupport positioning option');
            end
       end

    end
    
end






end