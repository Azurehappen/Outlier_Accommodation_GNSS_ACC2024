function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    LTSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,Option)
% xk is the prior state
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
yCov = p.sig_y^2.*eye(num);
E_R  = sqrt(yCov^(-1));
E_P  = sqrt(Pcov^(-1));
n    = size(H_os,2);
mu_x = zeros(n,1);
A    = [E_R*H_os; E_P];
c    = [E_R*res; E_P*mu_x];
h    = floor((num+n+n)/2)+floor((n+1)/2);
if nargin == 7
    switch Option
        case 1
            [rew,~] = ltsregres(A,c,'plots',0,'intercept',0);
        case 2             
             [rew,~] = ltsregres(A,c,'plots',0,'intercept',0,'h',h);
%             [rew,~] = ltsregres(A,c,'plots',0,'intercept',0);
    end 
else
    [rew,~] = ltsregres(A,c,'plots',0,'intercept',0);
end
        
%delta_x = rew.slope; 
b_LTS   = rew.flag; % binary vector corresponding to measurement selection
by      = b_LTS(1:num);
bx      = b_LTS(num+1:end);
nsv     = sum(by); % number of measurements used
dnsv    = num - nsv; % number of measurements discarded
nprior  = sum(bx); % number of priors used
dnprior = size(bx,1) - nprior; % number of priors discarded
% mslc_count = [nsv; dnsv; nprior; dnprior]; % msr selection count

Pby = diag(by);
Pbx = diag(bx);
PhiH = Pby*H_os;
delta_x = ((PhiH'*yCov^(-1)*PhiH+ Pbx'*Pcov^(-1)*Pbx)^(-1))*...
    (PhiH'*yCov^(-1)*Pby*res + Pbx'*Pcov^(-1)*Pbx*mu_x);
%-----------------------%
x_hat = xk + delta_x;
pos = x_hat(1:3);
clock_bias = x_hat(4);
err = norm(p.P_base - pos); % error check for debugging

% removes zero columns & rows to prevent singular matrix inversion
% when any sat constellation isn't present after measurement selection
Hbar    = PhiH;
zerocol = ~any(PhiH,1); % assigns 0 if col of 0s is present
Hbar(:,zerocol) = [];   % removes cols of 0s
zerorow = ~any(Hbar,2); % assigns 0 if row of 0s is present
Hbar(zerorow,:) = [];   % removes rows of 0s

% Compute GDOP
yCovbar = p.sig_y^2.*eye(size(Hbar,1));
GDOP   = sqrt(trace((Hbar'*yCovbar^(-1)*Hbar)^(-1)));

% alternative calculation to check if any different from delta_x
Pcov2 = diag(Pcov); Pcov2(zerocol) = []; Pcov2 = diag(Pcov2);
xk2 = xk; xk2(zerocol) = []; yCov2 = p.sig_y^2.*eye(size(Hbar,1));
mu_x(zerocol) = []; 
Pbx2 = diag(Pbx); Pbx2(zerocol) = []; Pbx2 = diag(Pbx2);
Pby_res = Pby*res; Pby_res(zerorow) = [];

delta_x2 = ((Hbar'*yCov2^(-1)*Hbar + Pbx2'*Pcov2^(-1)*Pbx2)^(-1))*...
    (Hbar'*yCov2^(-1)*Pby_res + Pbx2'*Pcov2^(-1)*Pbx2*mu_x);
%-----------------------%
x_hat2 = xk2 + delta_x2;
pos2 = x_hat2(1:3);
clock_bias2 = x_hat2(4);
err2 = norm(p.P_base - pos2);

% % g check
% g_LTS = rew.h;
% g1 = lts_g(A,num);
% 
% if abs(g_LTS-g1) ~= 1
%     diff = g_LTS-g1
%     fprintf('LTS: difference in g isn''t equal to 1 at epoch = %1.0f \n. ', p.idx)
% %     pause
% end

% compute cost
% A_b  = [E_R*Pby*H_os; E_P*Pbx];
% c_b  = [E_R*Pby*res; E_P*Pbx*mu_x];
% cost = (norm(A_b*delta_x - c_b))^2;
end