function [H_clk,x_clk] = formClkStatesAndH(num_sv)
% Create clock offset matrix for multi-GNSS
total = sum(num_sv);
num_sys = sum(num_sv~=0);
H_clk = zeros(total,num_sys);
col_ind = 0;
row_ind = 1;
for i = 1:length(num_sv)
    if num_sv(i) == 0
        continue;
    end
    H_clk(col_ind+1:col_ind+num_sv(i),row_ind)=1;
    col_ind = col_ind + num_sv(i);
    row_ind = row_ind + 1;
end
x_clk = zeros(num_sys,1);

end