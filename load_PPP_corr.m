function [p,obs] = load_PPP_corr(p,data_base,IGS_name,eph,obs,USTEC_folderpath,code_bia)

% USTEC & IGS data
if (p.post_mode==1 || p.post_mode==3)% If PPP, parse the iono correction
    if ~isempty(USTEC_folderpath) && p.L2enable ~= 1
        matname = USTEC_folderpath + "USTEC.mat";
        if exist(matname,'file')==2
            load(matname);
            p.USTEC = USTEC;
        else
            load_USTEC(p.t(1),p.t(end),USTEC_folderpath);
            USTEC = parser_USTEC(USTEC_folderpath);  
            p.USTEC = USTEC;
            save(matname,'USTEC')
        end
    end
    
    if ~isempty(IGS_name)
        matname = [IGS_name '.mat'];
        if exist(matname,'file')==2
            load(matname);
            p.IGS = IGS;
        else
            IGS = parser_CNES(IGS_name);
            p.IGS = IGS;
            save(matname,'IGS')
        end
    else
        fprintf('No IGS data or load error');
    end
    
end

if ~isempty(code_bia)
    p.code_bia = parser_bia(p,code_bia);
    p.bia_type = 1;
end

if p.L2enable == 1 && p.bia_type == 1
    if p.enableGPS == 1
        Factor = p.L2freq^2/(p.L1freq^2-p.L2freq^2);
        [n,r] = size(obs.GPS(1).data.P);
        obs.Iono_GPS = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.GPS.bia_C1C(i);
            bia_L2 = p.code_bia.GPS.bia_C2L(i);
            if bia_L2 ~= 0 && bia_L1 ~= 0
                for j = 1:r-400
                    data_L1 = obs.GPS(1).data.P(i,j:j+400);
                    data_L2 = obs.GPS(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_GPS(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
    if p.enableGAL == 1
        Factor = p.E5bfreq^2/(p.E1freq^2-p.E5bfreq^2);
        [n,r] = size(obs.GAL(1).data.P);
        obs.Iono_GAL = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.GAL.bia_C1C(i);
            bia_L2 = p.code_bia.GAL.bia_C7Q(i);
            if bia_L2 ~= 0 && bia_L1 ~= 0
                for j = 1:r-400
                    data_L1 = obs.GAL(1).data.P(i,j:j+400);
                    data_L2 = obs.GAL(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_GAL(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
    if p.enableBDS == 1
        Factor = p.B2afreq^2/(p.B1freq^2-p.B2afreq^2);
        [n,r] = size(obs.BDS(1).data.P);
        obs.Iono_BDS = zeros(n,r);
        for i = 1:n
            bia_L1 = p.code_bia.BDS.bia_C2I(i);
            bia_L2 = p.code_bia.BDS.bia_C7I(i);
            if bia_L2 ~= 0 && bia_L1 ~= 0
                for j = 1:r-400
                    data_L1 = obs.BDS(1).data.P(i,j:j+400);
                    data_L2 = obs.BDS(2).data.P(i,j:j+400);
                    if isempty(find(data_L2 == 0, 1)) && isempty(find(data_L1 == 0, 1))
                        obs.Iono_BDS(i,j) = Factor * ...
                            mean(data_L2 - p.c*bia_L2*1e-9 - data_L1 + p.c*bia_L1*1e-9);
                    end
                end
            end
        end
    end
end
p.eph_b = eph; p.obs_b = [];
% Base station data
if p.post_mode==2 && ~isempty(data_base)
%     matname = [data_base '_nav.mat'];
%     navname = [data_base '.nav'];
    obsname = [data_base '.obs'];
%     if exist(matname,'file')==2 % Check if the data already been parsed
%         load(matname);
%         p.eph_b = eph_b;
%     else
%         % Get ephemeris data (.nav file, RINEX verion 3.03)
%         if exist(matname,'file')==2
%             eph_b = parser_eph(p,navname);
%             save ([data_base,'_nav.mat'], 'eph_b');
%             p.eph_b = eph_b;
%         else
%             p.eph_b = eph;
%         end
%     end
    matname = [data_base '_obs.mat'];
    if exist(matname,'file')==2 % Check if the data already been parsed
        load(matname);
        p.obs_b = obs_b;
    else
        % Get observables data (.obs file, RINEX verion 3.03)
        obs_b = parser_obs(obsname);
        save ([data_base,'_obs.mat'], 'obs_b');
        p.obs_b = obs_b;
    end
end

end