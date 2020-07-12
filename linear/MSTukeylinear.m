function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    MSTukeylinear(p,xk,Pcov,H_offset,s_pos_ecef,y,bisqConst)
% xk is the prior state
% warning('off');
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
BisqConst = 4.685;
if nargin == 7
    BisqConst = bisqConst;
end
[delta_x,Wsquared] = Mestimator(c,A,BisqConst);
w = sqrt(diag(Wsquared)); % weight vector
W = diag(w); % weight matrix
% removes zero columns & rows to prevent singular matrix inversion
% when any sat constellation isn't present after measurement selection
Abar = W*A;
zerocol = ~any(Abar,1);   % find cols of 0s
Abar(:,zerocol) = [];     % remove cols of 0s
% zerorow = ~any(Abar,2); % assigns 0 if row of 0s is present
% Abar(zerorow,:) = [];   % removes rows of 0sAbar = W*A;

delta_xbar = inv(Abar'*Abar)*(Abar'*W*c);

%-----------------------%
x_hat = xk(1:4) + delta_xbar(1:4);
pos   = x_hat(1:3);
clock_bias = x_hat(4);
err1  = norm(p.P_base - pos);

% truncate Abar to keep only the portion related to measurements
Abar_trunc = Abar(1:num,:);
zerocol2 = ~any(Abar_trunc,1); % find cols of 0s
Abar_trunc(:,zerocol2) = [];    % remove cols of 0s
% compute GDOP
GDOP = sqrt(trace(inv(Abar_trunc'*Abar_trunc)));

wy = w(1:num);
wx = w(num+1:end);
nonzero_wy = find(wy);  % find non-zero weights
nonzero_wx = find(wx);  % find non-zero weights
nsv = numel(nonzero_wy);% number of measurements used
dnsv = num - nsv;       % number of measurements discarded
nprior  = size(nonzero_wx,1);   % number of priors used
dnprior = size(wx,1) - nprior; % number of priors discarded

% if err1 > 5
%     p.idx
%     pause
% end

% Compute Cost
% Ew = sqrt(weights);
% cost = (norm(Ew'*(A*delta_x - c)))^2;
end
%%%%=======================================================================
function [xEstimate,Weights] = Mestimator(y,H,bisqk)

    Tol = 1e-4;    % IRWLS iteration tolerance
    Totaliter = 0; % total iterations before convergence
    BisqK = 4.685; % Default Bisquare Normalizing Constant
    
    if nargin == 3
        BisqK = bisqk;
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
        [xEst,W]  = calcWLS(u,y,H,BisqK);
        diff      = norm(xEst - xEst_old);
        xEst_old  = xEst;
        iterCount = iterCount + 1; 
    end  
    xEstimate = xEst;
    Weights   = W;
    Totaliter = iterCount;
end
%--------------------------------------------------------------------------
function [xEst,W] = calcWLS(u,y,H,Bisqk)
    % calculates Weighted Least Squares
    m = size(H,1);
    w = zeros(m,1);
    for i = 1:1:m
        w(i) = calcWeight(u(i),Bisqk);  % get weights for u_i
    end
    W    = diag(w);    
    xEst = inv(H'*W*H)*(H'*W*y);        % Weighted LS
end
%--------------------------------------------------------------------------
function SigEst = calcSigEst(res)
% calculates sigma_hat, estimate of scale for Regression M-Estimator

    MAD     = median(abs(res-median(res))); % Max Absolute Deviation (MAD)
    factor  = 0.675;      % MAD normalization factor
    SigEst  = MAD/factor; % Normalized MAD    
end
%--------------------------------------------------------------------------
function weight = calcWeight(u,BisqK)
    % Caluclates weights for Weighted Least Squares
    % BisqK = bisquare tuning constant
    % u = res/sig 

    if abs(u) <= BisqK
        weight = (1 - (u/BisqK)^2)^2;
    else
        weight = 0;
    end
end
%--------------------------------------------------------------------------