function cpt = elevaz_check(p,cpt,u_pos)
% computing the azimuth and elevation angle from user to satellite (s)
% Remove the satellite values that not satisfy elevation require
%
% Input
%       p: paramters
%       cpt: computed data of svprn_mar k, corr_range, s_pos_ecef,num_sv
%       u_pos estimated user position
ind_prn = find(cpt.svprn_mark~=0);
num = length(cpt.num_sv);
count = 0;
for i = 1:num
    len = cpt.num_sv(i);
    for j = 1:len
        % Compute the elev and az
        [cpt.elev(count+j), cpt.az(count+j)] = sat_elev_azimuth(p,u_pos,cpt.sat_pos_Rcorr(:,count+j));
        if cpt.elev(count+j) < p.elev_mark
            % Eliminate the variables where not yield the elev mask,
            cpt.corr_range(count+j) = 0;
            cpt.svprn_mark(ind_prn(count+j)) = 0;
            cpt.prn_record(ind_prn(count+j)) = 0;
            cpt.num_sv(i) = cpt.num_sv(i)-1;
        end
    end
    count = count + len;
end
ind = find(cpt.corr_range==0);
if ~isempty(ind)
    cpt.corr_range(ind) = [];
    cpt.trop_delay(ind) = [];
    cpt.iono_delay(ind) = [];
    cpt.tp(ind)=[];
    cpt.s_pos_ecef(:,ind) = [];
    if p.post_mode == 1
        cpt.s_pos_prc(:,ind) = [];
    end
    cpt.s_v_ecef(:,ind) = [];
    cpt.sat_pos_Rcorr(:,ind) = [];
    cpt.sat_posprc_Rcorr(:,ind) = [];
    cpt.sat_v_Rcorr(:,ind) = [];
    cpt.elev(ind) = [];
    cpt.az(ind) = [];
end
end