function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    MSHuberlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,huberConst)
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
yCov  = p.sig_y^2.*eye(num);
E_R   = sqrt(yCov^(-1));
E_P   = sqrt(Pcov^(-1));
n     = size(H_os,2);
mu_x  = zeros(n,1);
A     = [E_R*H_os; E_P];
c     = [E_R*res; E_P*mu_x];
HuberConst = 1.345;
if nargin == 7
    HuberConst = huberConst;
end
[delta_x,Wsquared] = Mestimator(c,A,HuberConst);
w = sqrt(diag(Wsquared)); % vector of weights
%-----------------------%
x_hat = xk + delta_x;
pos   = x_hat(1:3);
clock_bias = x_hat(4);

% for GDOP computation, truncate A & weights to remove portion
% corresponding to prior & keep that corresponding to measurement.
A_trunc = A(1:num,:);                   % truncated A matrix
W_trunc = diag(w(1:num));   % weight matrix truncated
aug_A = W_trunc * A_trunc;  % 

% Remark: A_aug is likely to have 0's column after measurement selection by
% weights_trunc when all satellites observation of certain constellation
% are given 0 weight by Mestimator. This will lead to singular matrix 
% inversion in GDOP calculation. To prevent this, 0's column need to be 
% checked for and removed.

% remove 0's columns from A_aug (if present)
aug_Abar = aug_A;
zerocol  = ~any(aug_A,1); % index of 0's column
aug_Abar(:,zerocol) = []; % removes 0's column

% compute GDOP
GDOP = sqrt(trace(inv(p.sig_y^-2*(aug_Abar'*aug_Abar))));

wy = w(1:num);
wx = w(num+1:end);
nonzero_wy = find(wy); % find non-zero weights
nonzero_wx = find(wx); % find non-zero weights
nsv = numel(nonzero_wy); % number of measurements used
dnsv = num - nsv; % number of measurements discarded
nprior  = size(nonzero_wx,1); % number of priors used
dnprior = size(wx,1) - nprior; % number of priors discarded


% Compute Cost
% Ew = sqrt(weights);
% cost = (norm(Ew'*(A*delta_x - c)))^2;
end
%%%%=======================================================================
function [xEstimate,Weights] = Mestimator(y,H,const)

    Tol = 1e-4; % IRWLS iteration tolerance
    Totaliter = 0; % total iterations before convergence
    K = 1.345;     % tuning constant for Huber objective function
    
    if nargin == 3
        K = const;
    end    
    
    % Caclulate Ordinary LS
    OLS  = @(y,H) (H'*H)\H' * y; % Ordinary Least Squres
    xEst_old = OLS(y,H); % compute inital estimate
    diff = inf;   
    iterCount = 0;
    
    %Iteratively Reweighted Least Squares method 
    while diff >= Tol && iterCount < 1000       
        res    = y - H*xEst_old;
        sigEst = calcSigEst(res);
        u         = res./sigEst;        
        [xEst,W]  = calcWLS(u,y,H,K);
        diff      = norm(xEst - xEst_old);
        xEst_old  = xEst;
        iterCount = iterCount + 1; 
    end  
    xEstimate = xEst;
    Weights   = W;
    Totaliter = iterCount;
end
%--------------------------------------------------------------------------
function [xEst,W] = calcWLS(u,y,H,const)
    % const - normalizing constant for Huber objective function
    % calculates Weighted Least Squares
    m = size(H,1);
    w = zeros(m,1);
    for i = 1:1:m
        w(i) = calcWeight(u(i),const);  % get weights for u_i
    end
    W    = diag(w);    
    xEst = (H'*W*H)^(-1)*(H'*W*y);        % Weighted LS
end
%--------------------------------------------------------------------------
function SigEst = calcSigEst(res)
% calculates sigma_hat, estimate of scale for Regression M-Estimator

    MAD     = median(abs(res-median(res))); % Max Absolute Deviation (MAD)
    factor  = 0.675;      % MAD normalization factor
    SigEst  = MAD/factor; % Normalized MAD    
end
%--------------------------------------------------------------------------
function weight = calcWeight(u,const)
    % Caluclates weights for Weighted Least Squares
    % const - normalizing constant
    % u     - res/sig 

    if abs(u) <= const
        weight = 1;
    else
        weight = const/abs(u);
    end
end
%--------------------------------------------------------------------------