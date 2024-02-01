function [x_plus,cov_plus,J_ned,aug_cost] = mapUpdate(by, x_prior, cov_prior, y, H, R, Rot_e2g)
% Maximum A Posteriori state estimate and the augmented cost function
% when doing measurement selection
% INPUT: by  - measurement selection vector
%        E_P - sqrt(P^-1)
%        E_R - sqrt(R^-1)

E_P = chol(cov_prior^-1);    % cholesky decomp. of state prior info. matrix
E_R = chol(R^-1);    % cholesky decomp. of measurement info. matrix

if isempty(by)
    error('Input arg "by" is empty.');
else
    Phiby = diag(by);
    A = [E_R * Phiby * H; E_P];             
    c = [E_R * Phiby * y; zeros(length(x_prior),1)];   
    dx = (A'*A)^-1*A'*c; % LeastSquares % posterior state estimate
    aug_cost = norm(A*dx-c)^2;
end

x_plus = x_prior + Rot_e2g' * dx;
PhiH = Phiby * H;
J_ned    = PhiH' * R^-1 * PhiH + cov_prior^-1; % eqn. (13) in [1]
cov_plus = J_ned^-1;
cov_plus = Rot_e2g' * cov_plus * Rot_e2g;
end

