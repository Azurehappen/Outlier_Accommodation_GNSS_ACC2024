function [x, cov, flag] = checkClockReset(p, x0, cov0, num_user_states, res, cpt)

x = x0;
cov = cov0;
flag = false;
if sum(abs(res) > 5000) > 0.8*length(res)
    [estState,~,~] = weightLsSolver(p,cpt,true);
    [estState.vel,estState.clock_drift] = velSolver(estState.pos,cpt);
    x(num_user_states+1:end-1) = estState.clock_bias;
    x(end) = estState.clock_drift;
    cov(:, num_user_states+1:end) = zeros(length(x0),length(x0)-num_user_states);
    cov(num_user_states+1:end,:) = zeros(length(x0)-num_user_states,length(x0));
    clk_len = length(x0) - num_user_states - 1;
    cov(num_user_states+1:end-1, num_user_states+1:end-1) = 200^2*diag(ones(1,clk_len));
    cov(end, end) = 5^2;
    flag = true;
end

end