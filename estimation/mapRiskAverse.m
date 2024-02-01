function [flag,x_post,P_post,b,J_out,augcost,num_nodes,constraint] =...
    mapRiskAverse(num_constrain,y,H,P,R,J_l,x_prior,td_lambda,num_range)
% Solves RAPS using B&B integer optimization approach.
% Computes MAP state estimate using the selected measurements.
% OUTPUT:   x_post   - posterior state vector estimate
%           by       - measurement selection vector (binary)
%           augcost  - augmented cost for RAPS B&B
%           exitflag - see MATLAB function: intlinprog for description
% INPUT:    y - measurements
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

% see intlinprog function documentation for details
diagR  = diag(R);               % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;      % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
ieqLHS = -(sH.*sH)';            % inequality constraint LHS matrix (-G in the paper)
ieqLHS = ieqLHS(1:num_constrain,:);         % remove unconstrained rows
ieqRHS = diag(Jpminus - J_l);   % inequality constraint RHS vector (-g in the paper)
ieqRHS = ieqRHS(1:num_constrain,:);         % remove unconstrained rows
intcon = 1:1:m;                 % entries of b vector that can take only integer values
option = optimoptions(@intlinprog,'display','off'); % for output supression

constraint = -ieqLHS*ones(m,1)+diag(Jpminus(1:num_constrain,1:num_constrain));

E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

% ensure feasibility
b = ones(m,1);          % Initial b to be 1
num_nodes = 0;
% if ~all(ieqLHS*b <= ieqRHS)% && flag_rapid == true
%     loosen_count = 0;
%     while loosen_count <2 || ~all(ieqLHS*b <= ieqRHS)
%         max_neg_lhs = ieqLHS * b;
%         ind = find(ieqRHS < max_neg_lhs);
%         ieqRHS(ind) = 0.8 * max_neg_lhs(ind);  % loosen constraint
%         loosen_count = loosen_count + 1;
%     end
%     if ~all(ieqLHS*b <= ieqRHS)
%         warning('No constraint can be satisfied. Information Matrix Lower Bound too large');
%     end
% end
if ~all(ieqLHS*b <= ieqRHS) || m < n % m < 15?
    % Constraint cannot be satisfied, run TD.
    num_iter = 0;
    b = thresholdTest(td_lambda,P, y, H, R);
    [x_post,augcost] = MAP(b,y,H,E_R,E_P,x_prior);
    J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
    P_post = inv(J_out);
    flag = false;
    return;
end
flag = true;
xcurr = x_prior;
num_iter = 1;
total_trial = 15;
b = ones(m,1);
C = zeros(1, 2*total_trial);
C(1) = cost(xcurr,x_prior,P,b,H,y,R);
tic
while num_iter < total_trial
    Aeq = [];   % Coefficients for equality constrains: Aeq * b = Beq
    Beq = [];   % Values for equality constrains: Aeq * b = Beq
    res    = y-H*xcurr;             % residual
    s_res  = res./diagRs;           % residual scaled by meas. std (z-value)
    cost_b = (s_res).*(s_res);      % cost function coefficients for b optimization

    % solve for optimal integer measurement selection vector (Branch & Bound search)
    % solves the optimization problem in eqn (22) in [1].
    [b,~,exitflag,output] = intlinprog(cost_b,intcon,ieqLHS,ieqRHS,Aeq,Beq,lowerbound,upperbound,option);

    if (exitflag >= 1)
        % feasible solution is found
        num_nodes = num_nodes + output.numnodes;
        C(2*num_iter) = cost(xcurr,x_prior,P,b,H,y,R); % C(x_{l-1}, b_l)
        [x_post, augcost]  = MAP(b,y,H,E_R,E_P,x_prior); % compute x_l
    else
        % feasible solution not found 
        % Check the feasiblity of G*b before try.
        warning('feasible solution not found, this should not happen.')
    end
    
    xcurr = x_post;
    C(2*num_iter+1) = cost(xcurr,x_prior,P,b,H,y,R); % C(x_l, b_l)
    if abs(C(2*num_iter+1)-C(2*num_iter-1))<0.001
        break
    end
    num_iter = num_iter + 1;
end % while cnt_it

J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
%J_out % should meet spec
P_post = inv(J_out);

comp_t = toc; % computation time

% flag_pos = true;
% if flag == true
%     [flag_pos,flag_vel] = rapsValidation(td_lambda+1,P,y,H,R,b,num_range);
% end
% if flag_pos == false || flag_vel == false
%     num_iter = 0;
%     b = thresholdTest(td_lambda,P, y, H, R);
%     [x_post,augcost] = MAP(b,y,H,E_R,E_P,x_prior);
%     J_out   = calcJb(b,H,R,Jpminus); % posterior state information matrix
%     P_post = inv(J_out);
%     flag = false;
% end

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