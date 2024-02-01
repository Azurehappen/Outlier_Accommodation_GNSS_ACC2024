function [x_plus, cov_plus] = ekfUpdate(x_prior, cov_prior, meas_res, H, R)

% Innovation (or residual) covariance
S = H * cov_prior * H' + R;
% Near-optimal Kalman Gain
Kk = cov_prior * H' * S^(-1);
delta_x = Kk * meas_res;
x_plus = x_prior + delta_x;
I = eye(size(cov_prior));
cov_plus = (I - Kk*H)*cov_prior*(I - Kk*H)' + Kk*R*Kk';