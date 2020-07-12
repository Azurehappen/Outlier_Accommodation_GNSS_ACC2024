function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    LSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y)
% xk is the prior state
num = length(y); % total number of measurements
H = zeros(num,4);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
Pcov = diag(Pcov);
nsv     = num; % number of measurements used
dnsv    = 0; % number of measurements discarded
nprior  = size(Pcov,1); % number of priors used
dnprior = 0; % number of priors discarded
%-----------------------%      
for j=1:num
    R(j)=norm(s_pos_ecef(:,j)-xk(1:3));
    V= (xk(1:3)-s_pos_ecef(:,j))'/R(j);
    H(j,:)=[V 1];  
    r(j) = R(j)+sagnac(p,s_pos_ecef(:,j),xk);
    if ~isempty(H_offset)
        ind = find(H_offset(j,:)==1);
        if ~isempty(ind)
            off(j) = xk(4+ind);
        end
    end
end
res = y - r - xk(4)-off;
H_os = [H,H_offset];
%-----------------------%
yCov = p.sig_y^2.*eye(num); % noise covariance
delta_x = ((H_os'*yCov^(-1)*H_os+Pcov^(-1))^(-1))*(H_os'*yCov^(-1)*res);
GDOP = sqrt(trace((H_os'*yCov^(-1)*H_os)^(-1)));
%------------------------%
x_hat = xk + delta_x;
pos = x_hat(1:3);
clock_bias = x_hat(4);


% Compute Cost
% n = size(H_os,2);
% mu_x = zeros(n,1);
% E_R = sqrt(yCov^(-1));
% E_P = sqrt(Pcov^(-1));
% A_b = [E_R*H_os; E_P];
% c_b = [E_R*res; E_P*mu_x];
% cost = (norm(A_b*delta_x - c_b))^2; 
end