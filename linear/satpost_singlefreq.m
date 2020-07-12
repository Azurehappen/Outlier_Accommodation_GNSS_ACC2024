function [satlog] = satpost_singlefreq(p,eph,obs,ind,len_prn,sys_type)
% This function is to compute the satellite positions and correct 
% the psedoranges by sat clock bias
% Input: 
%       p --parameters
%       eph --ephemeris data
%       obs --observables
%       ind --index of observables data
%       sys_type --The system that be computed, 'GPS', 'GAL', 'GLO' and 'BDS'
% Outpu:
%       satlog.svprn_mark -- Mark the sat prn that be computed
%       satlog.s_pos_ecef -- Satellite position in ECEF frame
%       satlog.corr_range -- corrected pseudorange
%----------------------------%
% Initialize
satlog.svprn_mark = zeros(len_prn,1);
satlog.prn_record = zeros(len_prn,1);
satlog.s_pos_ecef = zeros(3,len_prn);
satlog.corr_range = zeros(len_prn,1);
obs_tr = obs.tr_sow(ind);
switch(sys_type)
    case 'GPS'
        obs_range = obs.GPS.P1(:,ind);
        Strength = obs.GPS.S1(:,ind);
        message_duration = p.gps.message_duration;
        eph_info = eph.gps;
    case 'GLO'
        obs_range = obs.GLO.P1(:,ind);
        Strength = obs.GLO.S1(:,ind);
        message_duration = p.glo.message_duration;
        obs_tr = time_shift(obs_tr - p.glo.lps_gps); % Correct time diff from GPS time to GLO time
        eph_info = eph.glo;
    case 'GAL'
        obs_range = obs.GAL.P1(:,ind);
        Strength = obs.GAL.S1(:,ind);
        message_duration = p.gal.message_duration;
        obs_tr = time_shift(obs_tr - p.gal.lps_gps); % Correct time diff from GPS time to GAL time
        eph_info = eph.gal;
    case 'BDS'
        obs_range = obs.BDS.P1(:,ind);
        Strength = obs.BDS.S1(:,ind);
        message_duration = p.bds.message_duration;
        obs_tr = time_shift(obs_tr - p.bds.lps_gps); % Correct time diff from GPS time to BDS time
        eph_info = eph.bds;
end

for j = 1:len_prn
    if (obs_range(j)~=0)&&(~isnan(obs_range(j)))&& j<=size(eph_info.t_oc,1)
        %----------------------%This part is to find the message in eph data row index
        tidx = find(eph_info.t_oc(j,:)~=0);
        dtr = limit_tgps(obs_tr-eph_info.t_oc(j,tidx)); % Limit time (in seconds) to 1-week
        if p.run_mode==1 && ~isempty(dtr(dtr>=-message_duration&dtr<=message_duration))
            % Post processing, find the index in eph corresponding to obs
            tidx = tidx(abs(dtr)==min(abs(dtr)));
            tidx = tidx(end);
        elseif p.run_mode==0 && ~isempty(dtr(dtr>=0&dtr<=message_duration))
            % Real time mode, cannot use the future eph message
            tidx = tidx(dtr == min(dtr(dtr>=0)));
            tidx = tidx(end);
        else
            tidx = [];
        end
        %----------------------%
        % Check the signal strength and sv health (health message o means ok)
        if ~isempty(tidx) && Strength(j)>=p.sig_strg && eph_info.SV_health(j,tidx)==0
            tp_prime = obs_range(j)/p.c;
            switch(sys_type)
                % compute ephemeris to satellite position and clock bias
                case 'GPS'
                    [sat_pos_ecef, ~, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,obs_tr+p.dt,tp_prime,'GPS');
                    if dt_sv ~= 0
                        satlog.svprn_mark(j) = 1;satlog.prn_record(j) = j;
                    end
                case 'GLO'
                    [sat_pos_ecef, ~, dt_sv] = geph2pos(p,eph_info,j,tidx,obs_tr,tp_prime);
                    if dt_sv ~= 0
                        satlog.svprn_mark(j) = 2;satlog.prn_record(j) = j;
                    end
                case 'GAL'
                    [sat_pos_ecef, ~, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,obs_tr,tp_prime,'GAL');
                    if dt_sv ~= 0
                        satlog.svprn_mark(j) = 3;satlog.prn_record(j) = j;
                    end
                case 'BDS'
                    [sat_pos_ecef, ~, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,obs_tr,tp_prime,'BDS'); 
                    if dt_sv ~= 0
                        satlog.svprn_mark(j) = 4;satlog.prn_record(j) = j;
                    end
            end
            if dt_sv ~= 0
                satlog.s_pos_ecef(:,j) = sat_pos_ecef;
                satlog.corr_range(j) = obs_range(j)+p.c*dt_sv;
            end
        end
    end
end
satlog.num_sv = sum(satlog.svprn_mark~=0);
end