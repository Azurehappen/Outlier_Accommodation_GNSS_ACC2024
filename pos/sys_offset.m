function [H_offset,x_offset] = sys_offset(num_sys)
% Create offset matrix for multi-GNSS
% Input:
%       num_offset: 4-1 vector, the number of satellite for each system
%       num_sys: 
% Outout:
%       offset: partial H matrix for offset paramters
total = sum(num_sys);
num_offset = sum(num_sys~=0)-1;
ind = find(num_sys~=0);
H_offset = zeros(total,num_offset);
start = 0;
if ~isempty(H_offset) % If only one system, don't need offset
    for i = 1:num_offset
        start = start + num_sys(ind(i));
        H_offset(start+1:start+num_sys(ind(i+1)),i)=1;
    end
end
x_offset = zeros(num_offset,1); % offset parts for state vector
end