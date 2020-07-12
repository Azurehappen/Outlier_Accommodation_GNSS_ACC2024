function p = load_PPP_corr(p,data_base,USTEC_folderpath,IGS_name,eph)

% USTEC & IGS data
if p.post_mode==1 || p.post_mode==3% If PPP, parse the iono correction
    if ~isempty(USTEC_folderpath)
        matname = [USTEC_folderpath '.mat'];
        if exist(matname,'file')==2
            load(matname);
            p.USTEC = USTEC;
        else
            USTEC = parser_USTEC(USTEC_folderpath);  
            p.USTEC = USTEC;
            save(matname,'USTEC')
        end
    else
        fprintf('No USTEC data or load error');
    end
    
    if ~isempty(IGS_name)
        matname = [IGS_name '.mat'];
        if exist(matname,'file')==2
            load(matname);
            p.IGS = IGS;
        else
            IGS = parser_IGS(IGS_name);
            p.IGS = IGS;
            save(matname,'IGS')
        end
    else
        fprintf('No IGS data or load error');
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