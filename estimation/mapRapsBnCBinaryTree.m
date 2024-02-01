function [x_star,P_post,b_star,cost_star,J_out,num_node] =...
    mapRapsBnCBinaryTree(y,H,P,R,J_l,x_prior)
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
[opt.x_star,opt.cost_star] = MAP(opt.b_star,y,H,E_R,E_P,x_prior);

res    = y-H*opt.x_star;             % residual
s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
coe_b = (s_res).*(s_res);      % cost function coefficients for b optimization
option = optimoptions(@linprog,'display','off'); % for output supression
% Initialize b by relaxation.
[b_frac,~,exitflag] = linprog(coe_b,-G,-g,[],[],lowerbound,upperbound,option);
b_frac(b_frac < 0.5) = 0;
b_frac(b_frac >= 0.5) = 1;
while ~all(G*b_frac >= g) || sum(b_frac) < n
    % Find the first 0 in b
    idx = find(b_frac == 0, 1);
    % If no 0s are found, break out of the loop
    if isempty(idx)
        break;
    end
    % Set the found 0 to 1
    b_frac(idx) = 1;
end
             
opt.b_star = b_frac;
[opt.x_star,opt.cost_star] = MAP(opt.b_star,y,H,E_R,E_P,x_prior);
opt.num_node = 0;
opt = binaryBranchAndCut(opt,opt.b_star,1,false);
x_star = opt.x_star;
b_star = opt.b_star;
cost_star = opt.cost_star;
num_node = opt.num_node;
J_out   = calcJb(opt.b_star,H,R,Jpminus); % posterior state information matrix
%J_out should meet spec
P_post = inv(J_out);

function opt = binaryBranchAndCut(opt,b,i_b,isRightBranch)
% DFS
if i_b > length(b)
    return;
end

% No computation for right branch since it is duplicate.
if isRightBranch
    % Go left
    b(i_b) = 0;
    opt = binaryBranchAndCut(opt,b,i_b+1,true);
    % Go right
    b(i_b) = 1;
    opt = binaryBranchAndCut(opt,b,i_b+1,false);
end

if sum(b) < n || ~all(G*b_frac >= g)
    % If index exceed or current node does not satisfy constraints.
    return;
end

opt.num_node = opt.num_node + 1;
[x_post,cost_new] = MAP(b,y,H,E_R,E_P,x_prior);
if cost_new < opt.cost_star
    % Update best solution
    opt.x_star = x_post;
    opt.b_star = b;
    opt.cost_star = cost_new;
end

% Go left
b(i_b) = 0;
opt = binaryBranchAndCut(opt,b,i_b+1,false);
% Go right
b(i_b) = 1;
opt = binaryBranchAndCut(opt,b,i_b+1,true);

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