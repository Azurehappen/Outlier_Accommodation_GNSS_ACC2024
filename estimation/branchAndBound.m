function opt = branchAndBound(opt,f,A,g,lowerbound,upperbound,option)

if ~all(A*upperbound <= g)
    return;
end

[b_frac,cost,exitflag] = linprog(f,A,g,[],[],lowerbound,upperbound,option);
opt.num_node = opt.num_node + 1;
if exitflag ~= 1 || cost >= opt.cost_star
    return;
end

is_all_integer = all(abs(b_frac) <= 1e-5 | abs(b_frac-1) <= 1e-5);
if is_all_integer == true
    b = b_frac;
    if cost < opt.cost_star
        % Replace values close to 0 with 0
        b(abs(b) <= 1e-5) = 0;
        % Replace values close to 1 with 1
        b(abs(b-1) <= 1e-5) = 1;
        opt.b_star = b;
        opt.cost_star = cost;
        opt.exitflag = 1;
    end
    return;
end

% Find all non-integer in b
non_i = find(abs(b_frac) > 1e-5 & abs(b_frac-1) > 1e-5);

for i=1:length(non_i)
    ubound = upperbound;
    ubound(non_i(i)) = 0;
    opt = branchAndBound(opt,f,A,g,lowerbound,ubound,option);
    lbound = lowerbound;
    lbound(non_i(i)) = 1;
    opt = branchAndBound(opt,f,A,g,lbound,upperbound,option);
end

end