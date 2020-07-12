function [svprn_mark,prn_record,s_pos_ecef,corr_range,num_sv] = satpost_corrpsedR_singlefreq_parfor(p,eph,obs,ind,len_prn,sys_type)
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
        tp_prime = obs_range(j)/p.c;
        t_sv = obs_tr-tp_prime;
        tidx = ephtidx(eph_info.t_oc{j},t_sv,eph_info.SV_health(j,:),message_duration);
        %----------------------%
        % Check the signal strength and sv health (health message o means ok)
        if ~isempty(tidx) && Strength(j)>=p.sig_strg
            switch(sys_type)
                % compute ephemeris to satellite position and clock bias
                case 'GPS'
                    [sat, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,t_sv,'GPS');
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(j) = 1;satlog.prn_record(j) = j;
                    end
                case 'GLO'
                    [sat, dt_sv] = geph2pos(p,eph_info,j,tidx,t_sv,'GLO');
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(j) = 2;satlog.prn_record(j) = j;
                    end
                case 'GAL'
                    [sat, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,t_sv,'GAL');
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(j) = 3;satlog.prn_record(j) = j;
                    end
                case 'BDS'
                    [sat, dt_sv] = eph2pos(p,eph_info,obs,j,tidx,t_sv,'BDS'); 
                    if ~isnan(sat.pos_ecef(1))
                        satlog.svprn_mark(j) = 4;satlog.prn_record(j) = j;
                    end
            end
            if ~isnan(sat.pos_ecef(1))
                satlog.s_pos_ecef(:,j) = sat.pos_ecef;
                satlog.corr_range(j) = obs_range(j)+p.c*dt_sv;
            end
        end
    end
end
satlog.num_sv = sum(satlog.svprn_mark~=0);

%----------------------%
svprn_mark = satlog.svprn_mark;
prn_record = satlog.prn_record;
s_pos_ecef = satlog.s_pos_ecef;
corr_range = satlog.corr_range;
num_sv     = satlog.num_sv;
%----------------------%
end