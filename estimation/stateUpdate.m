function [p,estState,res] = stateUpdate(p, cpt, dt)

%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;
x_minus = p.state0;
num_user_states = p.modeToNumUserStates(p.state_mode);
[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
if length(x_clk) + num_user_states + 1 ~= length(p.state0)
    error('current No. of sys does not match the previous epoch');
end
%------------------%
y = cpt.corr_range;
p.num_sats_window = [p.num_sats_window(2:length(p.num_sats_window)), length(y)];
yv = [];
if p.state_mode == p.pva_mode
    yv = cpt.doppler;
    H_clk = [H_clk; zeros(length(yv),size(H_clk, 2))];
end
num = length(y) + length(yv); % The number of measurement
H = zeros(num,num_user_states+length(x_clk)+1);
H(:,num_user_states+1:num_user_states+length(x_clk)) = H_clk;
if p.state_mode == p.pva_mode
    H(length(y)+1:end,end) = ones(length(yv),1);
    res_v = zeros(length(yv),1);
end
Range = zeros(length(y),1);
r = zeros(length(y),1);

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end

s_v_ecef = cpt.sat_v_Rcorr;
for j=1:length(y)
    Range(j)=norm(s_pos_ecef(:,j)-x_minus(1:3));
    los = (x_minus(1:3)-s_pos_ecef(:,j))'/Range(j)+...
        [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
    H(j,1:3)=los;
    range_r = norm(cpt.sat_pos_Rcorr(:,j)-x_minus(1:3));
    los_r = (x_minus(1:3)-cpt.sat_pos_Rcorr(:,j))'/range_r;
    if p.state_mode == p.pva_mode
        H(j+length(y),4:6) = los_r;
        res_v(j) = yv(j) - los_r*(x_minus(4:6) - s_v_ecef(:,j));
    end
    r(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),x_minus(1:3));
end
H_os = H;
R = constructMeasNoise(p, cpt, dt); %cpt.elev, cpt.svprn_mark
% measurement residual
res = y - r - H_os(1:length(y),4:end)*x_minus(4:end);
[x_minus, p.state_cov, flag] = checkClockReset(p, x_minus, p.state_cov, ...
    num_user_states, res, cpt); 
if flag == true
    res = y - r - H_os(1:length(y),4:end)*x_minus(4:end);
end
if p.state_mode == p.pva_mode
    res_v = res_v - x_minus(end);
    Rdop = 2*eye(length(yv));
    R = [R,zeros(length(y),length(yv));
        zeros(length(yv),length(y)),Rdop];
    res_all=[res;res_v];
else
    res_all=res;
end
if p.i == 3439
    i=1;
end
% y - f(x0) = H (x - x0);
zk = res_all + H_os * x_minus;
switch p.est_mode
    case p.ekf_est
        [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
    % case p.td_est
    %     [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_all, H_os, R);
    case p.map_est
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9);
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * p.state_cov * Rot_e2g;
        [x_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(ones(num,1),x_minus,cov_prior,res_all,H_os,R,Rot_e2g);
        p.num_meas_used = num;
    case p.td_est
        flag_rapid = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
        if p.state_mode == p.pva_mode && flag_rapid == true
            p.state_cov = zeros(size(p.state_cov));
            p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
            p.state_cov(end,end) = 20^2;
        end
        b = thresholdTest(p.td_lambda,p.state_cov, res_all, H_os, R);
        % inds = b==1;
        % H_td = H_os(inds,:);
        % res_td = res_all(inds);
        % R_td = R(inds,inds);
        % [x_plus, cov_plus] = ekfUpdate(x_minus, p.state_cov, res_td, H_td, R_td);
        lla = ecef2lla(x_minus(1:3)', 'WGS84');
        R_e2g=computeRotForEcefToNed(lla');
        R_pva = [R_e2g, zeros(3,6);
            zeros(3,3), R_e2g, zeros(3,3);
            zeros(3,6), R_e2g];
        Rot_e2g = [R_pva, zeros(9,length(x_minus)-9);
            zeros(length(x_minus)-9, 9), eye(length(x_minus)-9)];
        H_os = H_os * Rot_e2g';
        cov_prior = Rot_e2g' * p.state_cov * Rot_e2g;
        [x_plus,cov_plus,p.infor_ned,p.augcost] = ...
            mapUpdate(b, x_minus, cov_prior, res_all, H_os, R, Rot_e2g);
        p.num_meas_used = sum(b);
    case p.raps_ned_est
        % Solve in NED frame
        lla_deg = ecef2lla(x_minus(1:3)', 'WGS84');
        R_eg=computeRotForEcefToNed(lla_deg);
        if p.state_mode == p.pva_mode
            R_pva = [R_eg, zeros(3,6);
                zeros(3,3), R_eg, zeros(3,3);
                zeros(3,6), R_eg];
            flag_rapid = checkRapidNumSatChange(p.num_sats_window, sum(cpt.num_sv~=0));
            % Rapid num of sat detected, may entering an open sky area,
            % reset prior covariance
            if flag_rapid == true
                p.state_cov = zeros(size(p.state_cov));
                p.state_cov(1:end-1,1:end-1) = diag(150^2*ones(1,length(x_minus)-1));
                p.state_cov(end,end) = 20^2;
            end
        else
            R_pva = R_eg;
        end
        Rot_e2g = [R_pva, zeros(num_user_states,length(x_minus)-num_user_states);
            zeros(length(x_minus)-num_user_states, num_user_states), eye(length(x_minus)-num_user_states)];
        Ht = H_os * Rot_e2g';
        xt_minus = Rot_e2g*(x_minus - x_minus);
        Pt_minus = Rot_e2g*p.state_cov*Rot_e2g';
        if p.state_mode == p.pva_mode
            num_constrain = 6;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec;...
                p.raps.velhor_cov_spec; p.raps.velhor_cov_spec;...
                p.raps.velver_cov_spec]);
            p_clk = diag([p.raps.va_cov_spec*ones(3,1);...
                p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(6, length(x_minus)-6);
                zeros(length(x_minus)-6, 6), p_clk];
        elseif p.state_mode == p.pos_mode
            num_constrain = 3;
            cov_spec_ecef = diag([p.raps.poshor_cov_spec; ...
                p.raps.poshor_cov_spec; p.raps.posver_cov_spec]);
            p_clk = diag([p.raps.clk_cov_spec*ones(length(x_clk),1);...
                p.raps.dclk_cov_spec]);
            p_u = [cov_spec_ecef, zeros(3, length(x_minus)-3);
                zeros(length(x_minus)-3, 3), p_clk];
        end
        J_l = p_u^(-1);
        td_lambda_raps = 2;
        [flag,x_ned,cov_ned,b,J_out,p.augcost,num_nodes,constraint] = ...
            mapRiskAverse(num_constrain,res_all,Ht,Pt_minus,R,...
            diag(diag(J_l)),xt_minus,td_lambda_raps,length(y));
        cov_plus = Rot_e2g' * cov_ned * Rot_e2g;
        p.num_meas_used = sum(b);
        if p.state_mode == p.pva_mode
            p.raps_num_sat = sum(b(1:length(y)));
            % flag_pos = true;
            % if flag == true
            %     [flag_pos,flag_vel] = rapsValidation(p.td_lambda+1,Pt_minus, res_all, Ht, R, b, length(y));
            % end
            % % % Compute vel in NED
            % % prior_vel_ned = R_eg*x_minus(4:6);
            % % vel_ned = x_ned(4:6)+prior_vel_ned;
            % if flag_pos == false%... || flag_vel == false...
            %     %|| velocityValidation(prior_vel_ned,vel_ned,cov_ned(6,6)) == false
            %     cov_plus(1:6,1:6) = cov_plus(1:6,1:6) + diag(200^2*ones(6,1));
            % % elseif flag_vel == false ||...
            % %         velocityValidation(prior_vel_ned,vel_ned,cov_ned(6,6)) == false
            % %     cov_ned(4:6,4:6) = cov_ned(4:6,4:6) + diag(50^2*ones(3,1));
            % %     cov_ned(num_user_states+2,num_user_states+2) =...
            % %         cov_ned(num_user_states+2,num_user_states+2) + 50^2;
            % end
        else
            p.raps_num_sat = sum(b);
        end
        p.infor_ned = J_out;
        p.raps_num_nodes = num_nodes;
        p.constraint = constraint;
        p.raps_flag = flag;
        x_plus = Rot_e2g'*x_ned + x_minus;
    otherwise
        error('Incorrect state estimation mode configuration');
end

p.state0 = x_plus;
p.state_cov = cov_plus;

% HH = diag(b)*H_os;
% HH([1,5,7],:) = [];
% HH = HH(:, 1:4);
% hSqrtInv = (HH'*HH)^(-1);
% PDOP = sqrt(trace(hSqrtInv(1:3,1:3)));
% if PDOP > 5
%     estState.pos = [];
%     return;
% end

estState.pos = x_plus(1:3);
if p.state_mode == p.pva_mode
    estState.vel = x_plus(4:6);
end
estState.clock_bias = x_plus(4);
estState.clock_drift = x_plus(end);

clk_est = x_plus(4:end-1);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end


