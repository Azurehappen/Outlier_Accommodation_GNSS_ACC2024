function [t_out] = time_shift(t) 
%Shift the gps seconds interval in case of the difference < 0
t_out = t;
if (t < 0)
    t_out = t + 604800;

end