function cpt = earth_rotation_corr(p,cpt,clock_bias)
% Earth rotation correction.
%
%
% Parameters:
%   p:        parameters
%   clock_bias: estimated receiver clock bias
%
% Return values:
%   cpt: satllite computation data
%
% Reference:
%   IS-GPS-200
%
%
% Approved for use by University of California Riverside, Department of
% Electrical Engineering, Robotics and Control Lab.
%
% This program carries no warranty, not even the implied
% warranty of merchantability or fitness for a particular purpose.
%


% computations
%--------------------------------------------------------------%
cpt.sat_pos_Rcorr = zeros(3,length(cpt.tp));
cpt.sat_posprc_Rcorr = zeros(3,length(cpt.tp));
cpt.sat_v_Rcorr = zeros(3,length(cpt.tp));
for i = 1:length(cpt.tp)
%     theta = p.omge * (cpt.tp(i)-clock_bias - 5e-3);
    theta = p.omge * (cpt.tp(i)-clock_bias);
%     theta = p.omge * (cpt.tp(i));
    % rotation matrix
    R = [ cos(theta)  sin(theta)  0;   
         -sin(theta)  cos(theta)  0;
            0               0       1];
    cpt.sat_pos_Rcorr(:,i) = R*cpt.s_pos_ecef(:,i); 
    if p.post_mode == 1
        cpt.sat_posprc_Rcorr(:,i)=R*cpt.s_pos_prc(:,i);
    end
    cpt.sat_v_Rcorr(:,i) = R*cpt.s_v_ecef(:,i);
end


end