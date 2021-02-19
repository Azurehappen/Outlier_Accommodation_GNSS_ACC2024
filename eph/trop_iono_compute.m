function cpt = trop_iono_compute(p,eph,cpt,obs,re_pos,tdoy,ustec_i,user_t)
% Compute tropospheric delay and iono delay for the measurements

[p.lat, p.lon, p.h_r, ~, ~, ~] = ecef2llh(p,re_pos);
%%%%% convert geodetic height to orthometric height
%%%%% website: https://www.mathworks.com/matlabcentral/answers/97079-how-can-i-extract-the-orthometric-height-from-the-ellpsoiidal-height-in-the-mapping-toolbox-2-5-r20
load geoid;
%%%%% geodetic ellipsoidal separation
N = ltln2val(geoid, geoidrefvec, p.lat*180/pi, p.lon*180/pi);
%%%%% orthometric height of the receiver (Groves 2.121)
H_r=p.h_r-N;
% computing tropospheric delay for the reciever (s) using various
len = length(cpt.corr_range);
ind_prn = find(cpt.svprn_mark~=0);
for i = 1:len    
    % tropo delay (meter) computation using UNB3M model
%     [cpt.trop_delay(i), ~, ~, ~, cpt.IoFac(i)]=UNB3M(p.lat,H_r,tdoy,cpt.elev(i));
    
    % ecef2llh longitude range is [-pi,+pi] but IGGTrop takes longitude input
    % range is [0,2pi]. This mapping is done here.
    if p.lon <= 0
        longitude = -p.lon;
    else
        longitude = -p.lon + 2*pi;     
    end        
    % tropo delay (meter) computation using IGGTrop model
    % Reference Paper: IGGtrop_SH & IGGtrop_rH: Two Improved Empirical
    % Tropospheric Delay Models Based on Vertical Reduction Functions 
    IGGtrop_ZenithTropDelay = IGGtropSH_bl(rad2deg(longitude),rad2deg(p.lat),H_r/1000,tdoy);
    cpt.trop_delay(i) = (1.001/sqrt(0.002001 + sin(cpt.elev(i))^2))*IGGtrop_ZenithTropDelay;
    
    % Iono data from USTEC: https://www.ngdc.noaa.gov/stp/iono/ustec/products/    
    %---------------------------%
    % Select the frequncy
    switch cpt.svprn_mark(ind_prn(i))
        case 1 % GPS
            freq = p.L1freq;
%             switch obs.GPS.f1
%                 case 'L1'
%                     freq = p.L1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 2 % GLO
            freq = p.L1freq;
%             switch obs.GLO.f1
%                 case 'G1'
%                     freq = p.L1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 3 % GAL
            freq = p.E1freq;
%             switch obs.GAL.f1
%                 case 'E1'
%                     freq = p.E1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
        case 4 % BDS
            freq = p.B1freq;
%             switch obs.BDS.f1
%                 case 'B1'
%                     freq = p.B1freq;
%                 otherwise
%                     warning('Frequency type not support Iono delay computation, No result in this case');
%             end
    end
    %---------------------------%
    if ~isempty(ustec_i)
        % Computing Iono delay
        [cpt.iono_delay(i)] = ustec_iono_delay_computation(p,ustec_i,cpt.elev(i),cpt.az(i),user_t,freq);
    else
        if ~isempty(eph.ionoParameters)
            [cpt.iono_delay(i)] = klobuchar_model(p,eph.ionoParameters,cpt.elev(i),cpt.az(i),user_t.sow);
        else
            cpt.iono_delay(i) = 0;
        end
    end  
end