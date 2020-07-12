function  log = linear_gnss_ecef_parfor(p,eph,obs)
% This function is to implement GNSS positioning with
% linear mode (without Iono, Trop, Es correction)
% Output: log is a data struct that store the results.
%----------------------------------%
N = length(obs.tr_sow);         % The number of positioning points
p.MStkn = numel(p.MStkConst); 
p.MShbn = numel(p.MShbConst); 
p.LTSn  = numel(p.LTSOption);
p.TDn   = numel(p.TDLambda);
p.LSSn  = numel(p.LSSLambda);
pos_ecef   = cell(1,N); % Estimated position in ECEF
rover_clk  = cell(1,N); % Receiver clock bias (meter)
sv_num_GPS = NaN(1,N);  % The amount of GPS satellite be used
sv_num_GLO = NaN(1,N);  % The amount of GLO satellite be used
sv_num_GAL = NaN(1,N);  % The amount of GAL satellite be used
sv_num_BDS = NaN(1,N);  % The amount of BDS satellite be used
%--------------------------------------------------------------------------
if ~isempty(obs.GPS)
    num_obs_gps = size(obs.GPS.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_gps = 0;
end
if ~isempty(obs.GLO)
    num_obs_glo = size(obs.GLO.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_glo = 0;
end
if ~isempty(obs.GAL)
    num_obs_gal = size(obs.GAL.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_gal = 0;
end
if ~isempty(obs.BDS)
    num_obs_bds = size(obs.BDS.P1,1); % The maximum of PRN recorded in obs data
else
    num_obs_bds = 0;
end
res_GPS = NaN(num_obs_gps,N); % The residual at the end
res_GLO = NaN(num_obs_glo,N); % The residual at the end
res_GAL = NaN(num_obs_gal,N); % The residual at the end
res_BDS = NaN(num_obs_bds,N); % The residual at the end
max_num_sv = size([res_GPS; res_GAL; res_GLO;res_BDS],1);

% Initialize output
err_LS = NaN(1,N); % The position (Norm) error between estimated pos and true pos
hor_err_LS = NaN(1,N); % The horizontal position (Norm) error between estimated pos and true pos
ned_err_LS = NaN(1,N); % NED frame error
GDOPLS = NaN(1,N); % GDOP function (meter)
resLStemp = NaN(max_num_sv,1);
nsvLS = NaN(1,N); dnsvLS = NaN(1,N); npriorLS = NaN(1,N); dnpriorLS = NaN(1,N);

err_LSS = NaN(p.LSSn,N); hor_err_LSS = NaN(p.LSSn,N); 
ned_err_LSS = NaN(p.LSSn,N); GDOPLSS = NaN(p.LSSn,N);
resLSS = cell(p.LSSn,1); resLSStemp = NaN(max_num_sv*p.LSSn,1);
nsvLSS = NaN(p.LSSn,N); dnsvLSS = NaN(p.LSSn,N); npriorLSS = NaN(p.LSSn,N); dnpriorLSS = NaN(p.LSSn,N);

err_MShb = NaN(p.MShbn,N); hor_err_MShb = NaN(p.MShbn,N); 
ned_err_MShb = NaN(p.MShbn,N); GDOPMShb = NaN(p.MShbn,N);
resMShb = cell(p.MShbn,1); resMShbtemp = NaN(max_num_sv*p.MShbn,1);
nsvMShb = NaN(p.MShbn,N); dnsvMShb = NaN(p.MShbn,N); npriorMShb = NaN(p.MShbn,N); dnpriorMShb = NaN(p.MShbn,N);

err_MStk = NaN(p.MStkn,N); hor_err_MStk = NaN(p.MStkn,N); 
ned_err_MStk = NaN(p.MStkn,N); GDOPMStk = NaN(p.MStkn,N);
resMStk = cell(p.MStkn,1); resMStktemp = NaN(max_num_sv*p.MStkn,1);
nsvMStk = NaN(p.MStkn,N); dnsvMStk = NaN(p.MStkn,N); npriorMStk = NaN(p.MStkn,N); dnpriorMStk = NaN(p.MStkn,N);

err_LTS = NaN(p.LTSn,N); hor_err_LTS = NaN(p.LTSn,N);  
ned_err_LTS = NaN(p.LTSn,N); GDOPLTS = NaN(p.LTSn,N);
resLTS = cell(p.MShbn,1); resLTStemp = NaN(max_num_sv*p.LTSn,1);
nsvLTS = NaN(p.LTSn,N); dnsvLTS = NaN(p.LTSn,N); npriorLTS = NaN(p.LTSn,N); dnpriorLTS = NaN(p.LTSn,N);

err_TD = NaN(p.TDn,N); hor_err_TD = NaN(p.TDn,N); 
ned_err_TD = NaN(p.TDn,N); GDOPTD = NaN(p.TDn,N);
resTD = cell(p.TDn,1); resTDtemp = NaN(max_num_sv*p.TDn,1);
nsvTD = NaN(p.TDn,N); dnsvTD = NaN(p.TDn,N); npriorTD = NaN(p.TDn,N);dnpriorTD= NaN(p.TDn,N);

% x0 = p.x0;
for i = 1:N %N
    p.idx = i;
    % Mark the sat prn that be computed
    gpslog_svprn_mark = []; glolog_svprn_mark = [];
    gallog_svprn_mark = []; bdslog_svprn_mark = [];
    % Record the prn of each system for whose satllite been used
    gpslog_prn_record = []; glolog_prn_record = [];
    gallog_prn_record = []; bdslog_prn_record = [];
    % Satellite position in ECEF frame
    gpslog_s_pos_ecef = []; glolog_s_pos_ecef = [];
    gallog_s_pos_ecef = []; bdslog_s_pos_ecef = [];
    % corrected pseudorange
    gpslog_corr_range = []; glolog_corr_range = [];
    gallog_corr_range = []; bdslog_corr_range = [];
    % The number of satellite be computed
    gpslog_num_sv = 0; glolog_num_sv = 0;
    gallog_num_sv = 0; bdslog_num_sv = 0;

    cpt = struct();
    if ~isempty(obs.GPS)&& p.freq==1 && p.enableGPS == 1
        % GPS satellite position computation, Single frenquency receiver mode    
        [gpslog_svprn_mark,gpslog_prn_record,gpslog_s_pos_ecef,gpslog_corr_range,gpslog_num_sv] =...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_gps,'GPS');
    end
    if ~isempty(obs.GLO)&& p.freq==1 && p.enableGLO == 1
        % GLO satellite position computation, Single frenquency receiver mode
        [glolog_svprn_mark,glolog_prn_record,glolog_s_pos_ecef,glolog_corr_range,glolog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_glo,'GLO');
    end
    if ~isempty(obs.GAL)&& p.freq==1 && p.enableGAL == 1
        % GAL satellite position computation, Single frenquency receiver mode
        [gallog_svprn_mark,gallog_prn_record,gallog_s_pos_ecef,gallog_corr_range,gallog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_gal,'GAL');
    end
    if ~isempty(obs.BDS)&& p.freq==1 && p.enableBDS == 1
        % BDS satellite position computation, Single frenquency receiver mode
        [bdslog_svprn_mark,bdslog_prn_record,bdslog_s_pos_ecef,bdslog_corr_range,bdslog_num_sv] = ...
            satpost_corrpsedR_singlefreq_parfor(p,eph,obs,i,num_obs_bds,'BDS');
    end
    cpt.prn_record = [gpslog_prn_record;glolog_prn_record;gallog_prn_record;bdslog_prn_record];
    cpt.svprn_mark = [gpslog_svprn_mark;glolog_svprn_mark;gallog_svprn_mark;bdslog_svprn_mark];
    cpt.corr_range = [gpslog_corr_range;glolog_corr_range;gallog_corr_range;bdslog_corr_range];
    ind = find(cpt.corr_range==0);
    cpt.corr_range(ind) = [];
    cpt.s_pos_ecef = [gpslog_s_pos_ecef,glolog_s_pos_ecef,gallog_s_pos_ecef,bdslog_s_pos_ecef];
    cpt.s_pos_ecef(:,ind) = [];
    cpt.num_sv = [gpslog_num_sv,glolog_num_sv,gallog_num_sv,bdslog_num_sv];
    % elevation & azimuth
    cpt.elev = NaN(length(cpt.corr_range),1); cpt.az = NaN(length(cpt.corr_range),1);
    % trop delay and iono delay
    cpt.trop_delay = NaN(length(cpt.corr_range),1); cpt.iono_delay = NaN(length(cpt.corr_range),1);
    
    if sum(cpt.num_sv)>=p.min_sv
        index = round(p.Grdpos.t) == round(obs.tr_sow(i));
        grdpos = p.Grdpos.pos(:,index);
        if isempty(grdpos)
            continue;
        end

%         if p.state0(1)==0
%             [A,p.x0(4),~] = userpos(p,cpt.s_pos_ecef,cpt.corr_range,cpt.num_sv);
%         end
        % Check elevation
        p.x0(1:3) = grdpos;
        re_pos = p.x0(1:3);
        cpt = elevaz_check_linear(p,cpt,re_pos);
        if sum(cpt.num_sv)>=p.min_sv
            if ~isempty(p.eph_b) && ~isempty(p.obs_b)
                [cpt,n] = diff_corr_compute_linear(p,cpt,obs.tr_sow(i));
                if ~isempty(n) && sum(cpt.num_sv)>=p.min_sv
                    cpt = cpt_clear(cpt); % Clear the data where don't have diff correction
                    cpt.corr_range = cpt.corr_range - cpt.diff_corr;
                    %---------% Get prior for reciver clock bias
                    [~,x0_4,~] = userpos(p,cpt);                    
                    %---------%
                    [re_pos,clk_b,res,cpt.GDOP,cpt.nsv,cpt.dnsv,cpt.nprior,cpt.dnprior] = ...
                        userposlinear_parfor(p,cpt.s_pos_ecef,cpt.corr_range,cpt.num_sv,x0_4);                  
                    %--------------------save_result_linear---------------%
                    pos_ecef{i} = re_pos;
                    rover_clk{i} = clk_b;
                    sv_num_GPS(i) = cpt.num_sv(1);
                    sv_num_GLO(i) = cpt.num_sv(2);
                    sv_num_GAL(i) = cpt.num_sv(3);
                    sv_num_BDS(i) = cpt.num_sv(4);
                    ind_mark = cpt.svprn_mark ~= 0;
                    
                    % Find the obs_tr corresponding to the time in p.Grdpos.t 
                    % For example:
%                     index = round(p.Grdpos.t) == round(obs.tr_sow(i));
%                     grdpos = p.Grdpos.pos(:,index);
%                     if isempty(grdpos)
%                         continue;
%                     end

                    [ned_err_LS(:,i),hor_err_LS(:,i),err_LS(:,i),GDOPLS(:,i),...
                        nsvLS(i),dnsvLS(i),npriorLS(i),dnpriorLS(i),resLStemp(:,i)] = ...
                        save_errNorm_res_GDOP(1,grdpos,re_pos.LS,cpt.GDOP.LS,...
                        cpt.nsv.LS,cpt.dnsv.LS,cpt.nprior.LS,cpt.dnprior.LS,...
                        max_num_sv,ind_mark,res.LS);
                    
                    if p.eb_LSS == 1
                    [ned_err_LSS(:,i),hor_err_LSS(:,i),err_LSS(:,i),...
                        GDOPLSS(:,i),nsvLSS(:,i),dnsvLSS(:,i),npriorLSS(:,i),...
                        dnpriorLSS(:,i),resLSStemp(:,i)] =...
                        save_errNorm_res_GDOP(p.LSSn,grdpos,re_pos.LSS,...
                        cpt.GDOP.LSS,cpt.nsv.LSS,cpt.dnsv.LSS,cpt.nprior.LSS,...
                        cpt.dnprior.LSS,max_num_sv,ind_mark,res.LSS);
                    end
                    
                    if p.eb_MShb == 1
                    [ned_err_MShb(:,i),hor_err_MShb(:,i),err_MShb(:,i),...
                        GDOPMShb(:,i),nsvMShb(:,i),dnsvMShb(:,i),npriorMShb(:,i),...
                        dnpriorMShb(:,i),resMShbtemp(:,i)] = ...
                        save_errNorm_res_GDOP(p.MShbn,grdpos,re_pos.MShb,...
                        cpt.GDOP.MShb,cpt.nsv.MShb,cpt.dnsv.MShb,cpt.nprior.MShb,...
                        cpt.dnprior.MShb,max_num_sv,ind_mark,res.MShb);
                    end
                    
                    if p.eb_MStk == 1
                    [ned_err_MStk(:,i),hor_err_MStk(:,i),err_MStk(:,i),GDOPMStk(:,i),...
                        nsvMStk(:,i),dnsvMStk(:,i),npriorMStk(:,i),dnpriorMStk(:,i),...
                        resMStktemp(:,i)] = save_errNorm_res_GDOP(p.MStkn,...
                        grdpos,re_pos.MStk,cpt.GDOP.MStk,cpt.nsv.MStk,cpt.dnsv.MStk,...
                        cpt.nprior.MStk,cpt.dnprior.MStk,max_num_sv,ind_mark,res.MStk);
                    end
                    
                    if p.eb_LTS == 1
                    [ned_err_LTS(:,i),hor_err_LTS(:,i),err_LTS(:,i),GDOPLTS(:,i),...
                        nsvLTS(:,i),dnsvLTS(:,i),npriorLTS(:,i),dnpriorLTS(:,i),...
                        resLTStemp(:,i)] = save_errNorm_res_GDOP(p.LTSn,...
                        grdpos,re_pos.LTS,cpt.GDOP.LTS,cpt.nsv.LTS,cpt.dnsv.LTS,...
                        cpt.nprior.LTS,cpt.dnprior.LTS,max_num_sv,ind_mark,res.LTS);
                    end
                    
                    if p.eb_TD == 1
                    [ned_err_TD(:,i),hor_err_TD(:,i),err_TD(:,i),GDOPTD(:,i),...
                        nsvTD(:,i),dnsvTD(:,i),npriorTD(:,i),dnpriorTD(:,i),...
                        resTDtemp(:,i)] = save_errNorm_res_GDOP(p.TDn,...
                        grdpos,re_pos.TD,cpt.GDOP.TD,cpt.nsv.TD,cpt.dnsv.TD,...
                        cpt.nprior.TD,cpt.dnprior.TD,max_num_sv,ind_mark,res.TD);
                    end
                    %-----------------------------------------------------%                    
                end
            else
                warning('No differential source given')
            end
       end
    end    
end

% output
log.tr_prime  = obs.tr_prime; % utc time
log.gpst      = obs.tr_sow-obs.tr_sow(1);
log.pos_ecef  = pos_ecef;     % Estimated position in ECEF
log.rover_clk = rover_clk;    % Receiver clock bias (meter)

log.sv_num_GPS = sv_num_GPS; 
log.sv_num_GLO = sv_num_GLO; 
log.sv_num_GAL = sv_num_GAL; 
log.sv_num_BDS = sv_num_BDS;

log.err_LS  = err_LS; log.hor_err_LS = hor_err_LS;
log.ned_err_LS = ned_err_LS;
log.GDOPLS  = GDOPLS;
log.nsvLS = nsvLS;
log.dnsvLS = dnsvLS;
log.npriorLS = npriorLS;
log.dnpriorLS = dnpriorLS;
log.resLS = {resLStemp};

if p.eb_LSS == 1
    pointer1 = 1;
    for idx = 1:p.LSSn
        resLSS{idx} = resLSStemp(pointer1:pointer1+max_num_sv-1,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_LSS   = err_LSS; log.hor_err_LSS = hor_err_LSS;
    log.ned_err_LSS = ned_err_LSS;
    log.GDOPLSS = GDOPLSS;
    log.nsvLSS  = nsvLSS;
    log.dnsvLSS = dnsvLSS;
    log.npriorLSS = npriorLSS;
    log.dnpriorLSS = dnpriorLSS;
    log.resLSS = resLSS;
    log.LSSLambda = p.LSSLambda; log.LSSn = p.LSSn;
end

if p.eb_MShb == 1
    pointer1 = 1;
    for idx = 1:p.MShbn
        resMShb{idx} = resMShbtemp(pointer1:pointer1+max_num_sv-1,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_MShb  = err_MShb; log.hor_err_MShb = hor_err_MShb; 
    log.ned_err_MShb = ned_err_MShb;
    log.GDOPMShb = GDOPMShb;
    log.nsvMShb = nsvMShb;
    log.dnsvMShb = dnsvMShb;
    log.npriorMShb = npriorMShb;
    log.dnpriorMShb = dnpriorMShb;
    log.resMShb = resMShb;
    log.MShbConst = p.MShbConst; log.MShbn = p.MShbn;
end

if p.eb_MStk == 1
    pointer1 = 1;
    for idx = 1:p.MStkn
        resMStk{idx} = resMStktemp(pointer1:pointer1+max_num_sv-1,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_MStk  = err_MStk; log.hor_err_MStk = hor_err_MStk; 
    log.ned_err_MStk = ned_err_MStk;
    log.GDOPMStk = GDOPMStk;
    log.nsvMStk = nsvMStk;
    log.dnsvMStk = dnsvMStk;
    log.npriorMStk = npriorMStk;
    log.dnpriorMStk = dnpriorMStk;
    log.resMStk = resMStk;
    log.MStkConst = p.MStkConst; log.MStkn = p.MStkn;
end

if p.eb_LTS == 1
    pointer1 = 1;
    for idx = 1:p.LTSn
        resLTS{idx} = resLTStemp(pointer1:pointer1+max_num_sv-1,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_LTS   = err_LTS; log.hor_err_LTS = hor_err_LTS; 
    log.ned_err_LTS = ned_err_LTS;
    log.GDOPLTS = GDOPLTS;
    log.nsvLTS  = nsvLTS;
    log.dnsvLTS = dnsvLTS;
    log.npriorLTS = npriorLTS;
    log.dnpriorLTS = dnpriorLTS;
    log.resLTS = resLTS;
    log.LTSOption = p.LTSOption; log.LTSn = p.LTSn;
end

if p.eb_TD == 1
    pointer1 = 1;
    for idx = 1:p.TDn
        resTD{idx} = resTDtemp(pointer1:pointer1+max_num_sv-1,:);
        pointer1 = pointer1 + max_num_sv;
    end
    log.err_TD   = err_TD; log.hor_err_TD = hor_err_TD; 
    log.ned_err_TD = ned_err_TD; 
    log.GDOPTD = GDOPTD;
    log.nsvTD = nsvTD;
    log.dnsvTD = dnsvTD;
    log.npriorTD = npriorTD;
    log.dnpriorTD = dnpriorTD;
    log.resTD = resTD; 
    log.TDLambda = p.TDLambda; log.TDn = p.TDn;
end

end