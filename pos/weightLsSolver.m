function [estState,res,state_cov] = weightLsSolver(p,cpt,code_only_flag)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize
estState.clock_sys = dictionary;
estState.clock_sys(p.gps.sys_num) = NaN;
estState.clock_sys(p.glo.sys_num) = NaN;
estState.clock_sys(p.gal.sys_num) = NaN;
estState.clock_sys(p.bds.sys_num) = NaN;

x_pisb = [];
x_amb = [];
if p.post_mode == p.mode_rtkfloat && code_only_flag == false
    estState.phaseToGpsRho = dictionary;
    estState.phaseToGpsRho(p.gps.sys_num) = NaN;
    estState.phaseToGpsRho(p.glo.sys_num) = NaN;
    estState.phaseToGpsRho(p.gal.sys_num) = NaN;
    estState.phaseToGpsRho(p.bds.sys_num) = NaN;
    x_pisb = zeros(sum(cpt.num_sv ~= 0), 1);
    x_amb = zeros(sum(cpt.num_sv), 1);
end

x0 = p.state0(1:3);
[H_clk,x_clk] = formClkStatesAndH(cpt.num_sv);
xk = [x0;x_clk;x_pisb;x_amb];

y_rho = cpt.corr_range;
num = length(y_rho); % The number of measurement
y_phi = [];
if p.post_mode == p.mode_rtkfloat && code_only_flag == false
    y_phi = cpt.phase_m;
    if length(y_rho) ~= length(y_phi)
        warning('Num of pseudorange != phase');
    end
end
H = zeros(num+length(y_phi), length(xk));
H(1:num, 4:(4+length(x_clk)-1)) = H_clk;

H_pisb = zeros(length(y_phi),length(x_pisb));
if ~isempty(H_pisb)
    rowIdx = 1;  % To keep track of the current row in H
    colIdx = 1;  % To keep track of the current column in H
    for i = 1:length(cpt.num_sv)
        if cpt.num_sv(i) ~= 0
            H_pisb(rowIdx:rowIdx+cpt.num_sv(i)-1, colIdx) = 1;
            rowIdx = rowIdx + cpt.num_sv(i);
            colIdx = colIdx + 1;
        end
    end
    H(num+1:end, 5+length(x_isb):4+length(x_isb)+length(x_pisb)) = H_pisb;
    H(num+1:end, 5+length(x_isb)+length(x_pisb):end) = diag(cpt.wavelength);
end

s_pos_ecef = cpt.s_pos_ecef;
if p.post_mode == 1 && p.IGS_enable == 1
    s_pos_ecef = cpt.s_pos_prc;
end

Range = zeros(num,1);
geo_range = zeros(num+length(y_phi),1);
meas = [y_rho;y_phi];
if code_only_flag == true
    W = eye(num+length(y_phi));
else
    R = constructMeasNoise(p, cpt, 1);
    W = diag(1./diag(R));
end
for iter=1:p.Nls
    for j=1:num
        Range(j)=norm(s_pos_ecef(:,j)-xk(1:3));
        los = (xk(1:3)-s_pos_ecef(:,j))'/Range(j)+...
            [-s_pos_ecef(2,j)*p.omge/p.c s_pos_ecef(1,j)*p.omge/p.c 0];
        H(j,1:3)=los;
        geo_range(j) = Range(j)+sagnac(p,s_pos_ecef(:,j),xk);
        if p.post_mode == p.mode_rtkfloat && code_only_flag == false
            H(num+j,1:3)=los;
            geo_range(num+j) = geo_range(j);
        end
    end
    res = meas - geo_range - H(:, 4:end) * xk(4:end);
    delta_x = (H'*W*H)^(-1)*H'*W*(res);
    xk=xk+delta_x;
    if (norm(delta_x(1:4)) < p.LSthrsh)
        break;
    end
    if (iter>p.Nls)&& (norm(delta_x) > p.LSthrsh)
        warning('Postion path length iteration failed in user_pos calculation');
    end
end
%------------------%
estState.pos = xk(1:3);
estState.clock_bias = xk(4);
clk_est = xk(4:4+length(x_clk)-1);

state_cov = (H'*W*H)^(-1);
res = res(1:num);
j = 1;
for i = 1:length(cpt.num_sv)
    if cpt.num_sv(i) == 0
        continue;
    end
    estState.clock_sys(i) = clk_est(j);
    j=j+1;
end

end
