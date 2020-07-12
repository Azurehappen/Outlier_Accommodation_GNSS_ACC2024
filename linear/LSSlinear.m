function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    LSSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,lamda)
% xk is the prior state
% s is the estimated outliers
num = length(y); % total number of measurements
H = zeros(num,4);
R = zeros(num,1);
r = zeros(num,1);
off = zeros(num,1);
Pcov = diag(Pcov);
s = zeros(num,1);
s_old = ones(num,1);
nsv     = num; % number of measurements used
dnsv    = 0; % number of measurements discarded
nprior  = size(Pcov,1); % number of priors used
dnprior = 0; % number of priors discarded
while norm(s_old-s)>1e-3
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
H_os = [H,H_offset];
res = y - s - r - xk(4)-off;
%-----------------------%
yCov = p.sig_y^2.*eye(num); % noise covariance
GDOP = (H_os'*yCov^(-1)*H_os+Pcov^(-1))^(-1);
delta_x = GDOP*(H_os'*yCov^(-1)*res);
%-----------------------%
s_old = s;
Dr = y - r - xk(4)-off - H_os*delta_x;
s = outlier_est(Dr,diag(sqrt(yCov)),lamda);

end

GDOP = sqrt(trace((H_os'*yCov^(-1)*H_os)^(-1)));
%------------------------%
x_hat = xk + delta_x;
pos = x_hat(1:3);
clock_bias = x_hat(4);


function s = outlier_est(ri,sigma,lamda)
    s = zeros(length(ri),1);
    for i =1:length(ri)
        if ri(i)>sigma(i)*lamda
            s(i) = sigma(i)*(ri(i)-sigma(i)*lamda);
        elseif ri(i)<-sigma(i)*lamda
            s(i) = sigma(i)*(ri(i)+sigma(i)*lamda);
        else
            s(i) = 0;
        end
    end
end


end