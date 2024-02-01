function [x_star,P_post,b_star,cost_star,J_out,num_node] =...
    mapRapsBranchCut(y,H,P,R,J_l,x_prior,Rot_e2g)
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

diagR  = diag(R);           % diagonal entries of measurement covariance
diagRs = sqrt(diagR);
sH     = diag(diagRs) \ H;  % scale row i by {\sigma_i)^{-1}: diagRs^-1 * H
G = (sH.*sH)';              % G in the paper
G = G(1:3,:);               % remove unconstrained rows
g = diag(J_l - Jpminus(1:3,1:3));    % g in the paper
g = g(1:3,:);               % remove unconstrained rows
E_P = chol(Jpminus);    % cholesky decomp. of state prior info. matrix
E_R = chol(Jrminus);    % cholesky decomp. of measurement info. matrix

b_star = ones(m,1);          % Initial b to be 1
[x_star,cost_star,dx_star] = MAP(b_star,y,H,E_R,E_P,x_prior,Rot_e2g);
num_node = 0;
marker = []; % Marker to denote end of a BFS level
queue = {b_star};
queue{end+1} = marker;
visited = {mat2str(b_star)};

while ~isempty(queue)
    b = queue{1}; % Get the first element
    queue(1) = []; % Remove the first element

    % Check if we've hit the marker
    if isequal(b, marker)
        if isempty(queue)
            break;  % If there's no other node after marker, we're done
        end
        visited = {};  % Clear the visited list
        queue{end+1} = marker;  % Add marker to the end for next level
        continue;  % Skip the rest of the loop and process next node
    end

    % Branching: set one value in b to 0 and check the constraint
    for i = 1:m
        if b(i) == 0
            continue;
        end
        b_temp = b;
        b_temp(i) = 0;
        b_temp_str = mat2str(b_temp);

        if ismember(b_temp_str, visited)
            % Already visited this combination, so skip
            continue;
        end
        % Mark this combination as visited
        visited{end+1} = b_temp_str;

        if all(G*b_temp < g) || sum(b) < m-n
            % Constraint violated, so do not branch further from this vector
            continue;
        end
        
        num_node = num_node + 1;
        [x_post,cost_new,dx] = MAP(b_temp,y,H,E_R,E_P,x_prior,Rot_e2g);
        if cost_new < cost_star
            % Update best solution
            dx_star = dx;
            x_star = x_post;
            b_star = b_temp;
            cost_star = cost_new;
        end

         % Add to the queue for further branching
         queue{end+1} = b_temp;
    end
end
J_out   = calcJb(b_star,H,R,Jpminus); % posterior state information matrix
%J_out should meet spec
P_post = inv(J_out);

end

% function opt = binaryBranchAndBound(opt,coe_b,b,i_b,G,g,isLeftBranch)
% % DFS
% if i_b > length(b)
%     return;
% end
% % No computation for left branch since it is duplicate.
% if ~isLeftBranch
%     opt.node = opt.node + 1;
%     if G*b >= g
%         % Compute cost
%         cost_b = coe_b * b;
%         if cost_b < opt.cost_b
%             opt.cost_b = cost_b;
%             opt.b = b;
%         else
%             return;
%         end
%     end
% end
% % Go left
% b(i_b) = 0;
% opt = binaryBranchAndBound(opt,coe_b,b,i_b+1,G,g,true);
% % Go right
% b(i_b) = 1;
% opt = binaryBranchAndBound(opt,coe_b,b,i_b+1,G,g,false);
% 
% end

%--------------------------------------------------------------------------
function [x_post,aug_cost,dx] = MAP(by,y,H,E_R,E_P,x_prior, Rot_e2g)
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
    c = [E_R * Phiby * y; zeros(length(x_prior),1)];   % eqn. (10) in [1]
    dx = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*dx-c)^2;          % eqn. (11) in [1]
end
x_post = x_prior + Rot_e2g' * dx;

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