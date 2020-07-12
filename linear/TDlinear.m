function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,lambda)
% xk is the prior state
% lambda is threshold value
num = length(y); % total number of measurements
H = zeros(num,4);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
Pcov = diag(Pcov);
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
res = y- r - xk(4)-off;
H_os = [H,H_offset];
%-----------------------%
n = size(H_os,2);
mu_x = zeros(n,1);
yCov = p.sig_y^2.*eye(num);
reshat = H_os*mu_x;
resdiff = abs(res - reshat);

resCov = H_os*Pcov*H_os' + yCov; % Section VI:A, 1st sentence.
sig_ri = sqrt(diag(resCov));
ratio  = resdiff./sig_ri;
Lvec   = lambda.*ones(num,1); % vector of thresholds
by     = zeros(num,1); % binary decision vector

for i = 1:1:num
    % generate selection vector b_TD
    if ratio(i) <= Lvec(i) 
        by(i) = 1;     % keep res(i)
    else
        by(i) = 0;     % flag res(i)
    end
end
nsv = sum(by); % number of measurements used
dnsv    = num - nsv; % number of measurements discarded
nprior  = n; % number of priors used
dnprior = 0; % number of priors discarded

Pbx = eye(n);
Pby = diag(by);
PhiH  = Pby*H_os;
delta_x = ((PhiH'*yCov^(-1)*PhiH + Pbx'*Pcov^(-1)*Pbx)^(-1))*...
    (PhiH'*yCov^(-1)*Pby*res + Pbx'*Pcov^(-1)*Pbx*mu_x);
%-----------------------%
x_hat = xk + delta_x;
pos = x_hat(1:3);
clock_bias = x_hat(4);
err = norm(p.P_base - pos);

% removes zero columns & rows to prevent singular matrix inversion
% when any sat constellation isn't present after measurement selection
Hbar    = PhiH;
zerocol = ~any(PhiH,1);   % assigns 0 if col of 0s is present
Hbar(:,zerocol) = [];  % removes cols of 0s
zerorow = ~any(Hbar,2);   % assigns 0 if row of 0s is present
Hbar(zerorow,:) = [];  % removes rows of 0s

% Compute GDOP
yCovbar = p.sig_y^2.*eye(size(Hbar,1));
GDOP   = sqrt(trace((Hbar'*yCovbar^(-1)*Hbar)^(-1)));

% % alternative calculation to check if any different from delta_x
% Pcov2 = diag(Pcov); Pcov2(zerocol) = []; Pcov2 = diag(Pcov2);
% xk2 = xk; xk2(zerocol) = []; yCov2 = p.sig_y^2.*eye(size(Hbar,1));
% mu_x(zerocol) = []; Pbx = eye(sum(~zerocol));
% Pby_res = Pby*res; Pby_res(zerorow) = [];
% 
% delta_x2 = ((Hbar'*yCov2^(-1)*Hbar + Pbx'*Pcov2^(-1)*Pbx)^(-1))*...
%     (Hbar'*yCov2^(-1)*Pby_res + Pbx'*Pcov2^(-1)*Pbx*mu_x);
% %-----------------------%
% x_hat2 = xk2 + delta_x2;
% pos2 = x_hat2(1:3);
% clock_bias2 = x_hat2(4);
% err2 = norm(p.P_base - pos2);

% Compute Cost
% E_R  = sqrt(inv(yCov));
% E_P  = sqrt(inv(Pcov));
% A_b  = [E_R*Pby*H_os; E_P*Pbx];
% c_b  = [E_R*Pby*res; E_P*Pbx*mu_x];
% cost = (norm(A_b*delta_x - c_b))^2; 
end