function [pos,clock_bias,res] = userpos(p,cpt)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Measurement selection applied
% Input: 
%       s_pos_ecef: 3-by-N Satellite position in ECEF frame.
%       x0 : 3-by-1 initial interative coordinates in ECEF frame.
%       y: m-by-1 Corrected pseudorange.
%
% Output:
%       
%       
%       
%       

%-------------------%
% Initialize

x0 = p.state0;
[H_offset,x_offset] = sys_offset(cpt.num_sv);
xk = [x0;x_offset];
%------------------%
[pos,clock_bias,res] = LSsolver(p,xk,H_offset,cpt);
% switch p.select
%     case 0
%         [pos,clock_bias,res] = LSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 1
%         [pos,clock_bias,res,cost] = TDsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 2
%         [pos,clock_bias,res,cost] = LSSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 3
%         [pos,clock_bias,res,cost] = MSsolver(p,xk,H_offset,s_pos_ecef,y);
%     case 4
%         [pos,clock_bias,res,cost] = LTSsolver(p,xk,H_offset,s_pos_ecef,y);
% end

end
