function [x_star,P_post,b_star,cost_star,J_out,num_node] =...
    mapRapsBnBHeuristic(y,H,P,R,J_l,x_prior)
% Solves RAPS using B&B integer optimization approach.
% Computes MAP state estimate using the selected measurements.
% OUTPUT:   x_post   - posterior state vector estimate
%           by       - measurement selection vector (binary)
%           augcost  - augmented cost for RAPS B&B
%           exitflag - see MATLAB function: intlinprog for description
% INPUT:    res - measurements residuals: y - f(x0) = H*dx, where dx =
% x_plus - x_minus
%           H - measurement matrix
%           P - Prior Covariance matrix
%           r - Measurement Covariance matrix
%           J_l - Information Matrix Lower Bound
%           x_prior - prior state vector estimate
% reference: [1] - PPP_RAPS_Linear.pdf

Jpminus = P^-1;         % state prior info. matrix
Jrminus = R^-1;         % measurement info. matrix

[m,n] = size(H);
lowerbound = zeros(m,1); % lower bound on measurement selection vector b
upperbound = ones(m,1);  % upper bound on measurement selection vector b

diagR  = diag(R);           % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;  % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
G = (sH.*sH)';              % G in the paper
G = G(1:3,:);               % remove unconstrained rows
g = diag(J_l - Jpminus);    % g in the paper
g = g(1:3,:);               % remove unconstrained rows
E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

opt.b_star = ones(m,1);          % Initial b to be 1
if ~all(G*opt.b_star >= g)
    loosen_count = 0;
    while loosen_count <4 || ~all(G*opt.b_star >= g)
        min_Gb = G * opt.b_star;
        ind = find(g > min_Gb);
        g(ind) = 0.8 * min_Gb(ind);  % loosen constraint
        loosen_count = loosen_count + 1;
    end
end
if ~all(G*opt.b_star >= g) || m < n
    warning('No constraint can be satisfied. Information Matrix Lower Bound too large');
    b_star = opt.b_star;
    num_node = 0;
    [x_star,cost_star] = MAP(opt.b_star,y,H,E_R,E_P,x_prior);
    J_out   = calcJb(opt.b_star,H,R,Jpminus); % posterior state information matrix
    %J_out should meet spec
    P_post = inv(J_out);
    return;
end

num_node = 0;
xcurr = x_prior;
num_iter = 1;
total_trial = 20;
b = ones(m,1);
C = zeros(1, 2*total_trial);
C(1) = cost(xcurr,x_prior,P,b,H,y,R);
option = optimoptions(@linprog,'display','off'); % for output supression
opt = struct;
while num_iter < total_trial
    res    = y-H*xcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    coe_b = (s_res).*(s_res);      % cost function coefficients for b optimization
    opt.cost_star = Inf; %coe_b' * ones(m,1);
    opt.num_node = 0;
    opt.exitflag = 0;
    opt = branchAndBound(opt,coe_b,-G,-g,lowerbound,upperbound,option);

    if (opt.exitflag >= 1)
        % feasible solution is found
        num_node = num_node + opt.num_node;
        C(2*num_iter) = cost(xcurr,x_prior,P,opt.b_star,H,y,R); % C(x_{l-1}, b_l)
        %opt.cost_star = C(2*num_iter);
        [x_post, augcost]  = MAP(opt.b_star,y,H,E_R,E_P,x_prior); % compute x_l
    else
        % feasible solution not found
        % Check the feasiblity of G*b before try.
        warning('feasible solution not found, this should not happen.')
    end

    xcurr = x_post;
    C(2*num_iter+1) = cost(xcurr,x_prior,P,b,H,y,R);
    %opt.cost_star = C(2*num_iter+1); % C(x_l, b_l)
    if abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001
        break
    end
    num_iter = num_iter + 1;
end


J_out   = calcJb(opt.b_star,H,R,Jpminus); % posterior state information matrix
%J_out should meet spec
P_post = inv(J_out);

end

function [C] = cost(x,prior,P,b,H,y,R)
ex = x-prior;
Pi = inv(P);
ey = y - H*x;
Ri = inv(R);
Pb = diag(b);
if isempty(b)
    C1 = ex' * Pi * ex;
    C2 =  ey' *  Ri *  ey;
    C  = C1 + C2;
else
    C = ex' * Pi * ex + ey' * Pb * Ri * Pb * ey;
end
end

%--------------------------------------------------------------------------
function [x_post,aug_cost] = MAP(by,y,H,E_R,E_P,x_prior)
% Maximum A Posteriori state estimate and the augmented cost function
% when doing measurement selection
% INPUT: by  - measurement selection vector
%        E_P - sqrt(P^-1)
%        E_R - sqrt(R^-1)

if isempty(by)
    error('Input arg "by" is empty.');
else
    Phiby = diag(by);
    A = [E_R * Phiby * H; E_P];             % eqn. (10) in [1]
    c = [E_R * Phiby * y; E_P * x_prior];   % eqn. (10) in [1]
    x_post = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*x_post-c)^2;          % eqn. (11) in [1]
end
end
%--------------------------------------------------------------------------
function J = calcJb(by,H,R,Jpminus)
% Calculates posterior state information matrix when doing measurement selection
% INPUT: by - measurement selection vector
%        Jpminus - State information matrix prior, Jpminus = Pminus^-1

Pby  = diag(by); % Phi(by)
PhiH = Pby * H;
J    = PhiH' * R^-1 * PhiH + Jpminus; % eqn. (13) in [1]
end