function R = constructMeasNoise(p, cpt, dt)
% DGNSS measurement noise model from RTKLIB (rtkpos.c -> varerr)
% 2.0*(a*a+b*b/sinel/sinel+c*c)+d*d;
% a and b are both = code/carrier ratio (default: 300) * carrier-phase error STD (default: 0.003m)
% c is base length STD whose default is 0 for < 10km.
% d = speed of light * click stability (default: 5E-12) * delta time

elev_rad = cpt.elev;
sys_mark = cpt.svprn_mark;

R = zeros(2*length(elev_rad), 2*length(elev_rad));
fact_a = 300;
fact_b = 400;
fact_a_glo = 600;
fact_b_glo = 800;
sigma_a = 0.003;
sigma_b = 0.003;
code_a = fact_a*sigma_a;
code_b = fact_b*sigma_b;  
code_a_glo = fact_a_glo*sigma_a;
code_b_glo = fact_b_glo*sigma_b;
d = p.c * 5e-12 * dt;

for i = 1:length(elev_rad)
    if sys_mark(i) == p.glo.sys_num
        R(i,i) = measNoiseVar(p,elev_rad(i),cpt,i,code_a_glo,code_b_glo,d);
        R(length(elev_rad)+i,length(elev_rad)+i) =...
            measNoiseVar(p,elev_rad(i),cpt,i,sigma_a,sigma_b,d);
        continue;
    end
    R(i,i) = measNoiseVar(p,elev_rad(i),cpt,i,code_a,code_b,d);
    R(length(elev_rad)+i,length(elev_rad)+i) =...
            measNoiseVar(p,elev_rad(i),cpt,i,sigma_a,sigma_b,d);
end
if p.post_mode ~= p.mode_rtkfloat && p.post_mode ~= p.mode_rtkfix
    R = R(1:length(elev_rad),1:length(elev_rad));
end

    function var = measNoiseVar(p,elev_rad,cpt,i,sigma_a,sigma_b,d)
        if p.post_mode == p.mode_ppp
            var = p.ppp.sigma_cme^2+(cpt.iono_map_m(i)*p.ppp.sigma_iono)^2+sigma_a^2+d^2;
            %sigma_b^2/(sin(elev_rad))^2
        elseif p.post_mode == p.mode_dgnss
            var = 2*(sigma_a^2)+sigma_b^2 / (sin(elev_rad))^2 + d^2;
        end
    end

end