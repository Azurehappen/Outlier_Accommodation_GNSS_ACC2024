function [pos,clock_bias,res,GDOP] = userposlinear(p,s_pos_ecef,y,num_sys)
% This is solver for computing user position, receiver clock
% bias, satellite system offset.
% Linear mode
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
x0 = p.x0;
[H_offset,~] = sys_offset(num_sys);
ind = find(num_sys(2:4)~=0);
ISB = [p.ISBglo;p.ISBgal;p.ISBbds];
ISB_cov = [p.ISBglo_cov;p.ISBgal_cov;p.ISBbds_cov];
ISB = ISB(ind);ISB_cov = ISB_cov(ind);
xk = [x0;ISB]; % Prior
Pcov = [p.priorposcov;ISB_cov]; % covariance of prior state
xk(1:3) = xk(1:3)+sqrt(Pcov(1:3)).*randn(3,1);
%xk = xk + sqrt(Pcov) .* randn(length(Pcov),1);
%------------------%
[pos.LS,clock_bias.LS,res.LS,GDOP.LS] = LSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
if p.eb_LSS == 1
[pos.LSS,clock_bias.LSS,res.LSS,GDOP.LSS] = LSSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
end
if p.eb_MShb == 1
[pos.MShb,clock_bias.MShb,res.MShb,GDOP.MShb] = MSHuberlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
end
if p.eb_MStk == 1
[pos.MStk,clock_bias.MStk,res.MStk,GDOP.MStk] = MSTukeylinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
end
if p.eb_LTS == 1
[pos.LTS,clock_bias.LTS,res.LTS,GDOP.LTS] = LTSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
end

% p.lambda = 1; % threshold value for TD
% [pos.TD1,clock_bias.TD1,res.TD1,GDOP.TD1] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% p.lambda = 2; % threshold value for TD
% [pos.TD2,clock_bias.TD2,res.TD2,GDOP.TD2] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% p.lambda = 3; % threshold value for TD
% [pos.TD3,clock_bias.TD3,res.TD3,GDOP.TD3] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);

% switch p.select
%     case 0
%         [pos,clock_bias,res,cost] = LSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
%     case 1
%         [pos,clock_bias,res,cost] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
%     case 2
%         [pos,clock_bias,res,cost] = LSSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
%     case 3
%         [pos,clock_bias,res,cost] = MSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
%     case 4
%         [pos,clock_bias,res,cost] = LTSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% end
% 
end
