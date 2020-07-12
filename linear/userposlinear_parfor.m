function [pos,clock_bias,res,GDOP,nsv,dnsv,nprior,dnprior] = ...
    userposlinear_parfor(p,s_pos_ecef,y,num_sys,x0_4)
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
% if floor(p.idx) == 997
%     pause
% end
%-------------------%
% Initialize
x0 = p.x0;
x0(4) = x0_4;
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
[pos.LS{1},clock_bias.LS,res.LS{1},GDOP.LS,nsv.LS,dnsv.LS,nprior.LS,dnprior.LS]= ...
    LSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);

if p.eb_LSS == 1
pos.LSS  = cell(p.LSSn,1);  res.LSS = cell(p.LSSn,1);
GDOP.LSS = zeros(p.LSSn,1); clock_bias.LSS = zeros(p.LSSn,1);
nsv.LSS = zeros(p.LSSn,1); dnsv.LSS = zeros(p.LSSn,1);
nprior.LSS = zeros(p.LSSn,1);dnprior.LSS = zeros(p.LSSn,1);
for idx = 1:p.LSSn
    [pos.LSS{idx},clock_bias.LSS(idx),res.LSS{idx},GDOP.LSS(idx),...
    nsv.LSS(idx),dnsv.LSS(idx),nprior.LSS(idx),dnprior.LSS(idx)]= ...
    LSSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,p.LSSLambda(idx));
end
end

if p.eb_MShb == 1
pos.MShb  = cell(p.MShbn,1);  res.MShb = cell(p.MShbn,1);
GDOP.MShb = zeros(p.MShbn,1); clock_bias.MShb = zeros(p.MShbn,1);
nsv.MShb = zeros(p.MShbn,1); dnsv.MShb = zeros(p.MShbn,1);
nprior.MShb = zeros(p.MShbn,1);dnprior.MShb = zeros(p.MShbn,1);
for idx = 1:p.MShbn
    [pos.MShb{idx},clock_bias.MShb(idx),res.MShb{idx},GDOP.MShb(idx),...
    nsv.MShb(idx),dnsv.MShb(idx),nprior.MShb(idx),dnprior.MShb(idx)]=...
    MSHuberlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,p.MShbConst(idx));
end
end

if p.eb_MStk == 1
pos.MStk  = cell(p.MStkn,1);  res.MStk = cell(p.MStkn,1);
GDOP.MStk = zeros(p.MStkn,1); clock_bias.MStk = zeros(p.MStkn,1);
nsv.MStk = zeros(p.MStkn,1);dnsv.MStk = zeros(p.MStkn,1);
nprior.MStk = zeros(p.MStkn,1);dnprior.MStk = zeros(p.MStkn,1);
for idx = 1:p.MStkn
    [pos.MStk{idx},clock_bias.MStk(idx),res.MStk{idx},GDOP.MStk(idx),...
    nsv.MStk(idx),dnsv.MStk(idx),nprior.MStk(idx),dnprior.MStk(idx)]= ...
    MSTukeylinear(p,xk,Pcov,H_offset,s_pos_ecef,y,p.MStkConst(idx));
end
end


if p.eb_LTS == 1
pos.LTS  = cell(p.LTSn,1);  res.LTS = cell(p.LTSn,1);
GDOP.LTS = zeros(p.LTSn,1); clock_bias.LTS = zeros(p.LTSn,1);
nsv.LTS = zeros(p.LTSn,1);dnsv.LTS = zeros(p.LTSn,1);
nprior.LTS = zeros(p.LTSn,1);dnprior.LTS = zeros(p.LTSn,1);
for idx = 1:p.LTSn
    [pos.LTS{idx},clock_bias.LTS(idx),res.LTS{idx},GDOP.LTS(idx),...
    nsv.LTS(idx),dnsv.LTS(idx),nprior.LTS(idx),nprior.LTS(idx)] =...
    LTSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,p.LTSOption(idx)); 
end
end

if p.eb_TD == 1
pos.TD  = cell(p.TDn,1);  res.TD = cell(p.TDn,1);
GDOP.TD = zeros(p.TDn,1); clock_bias.TD = zeros(p.TDn,1);
nsv.TD = zeros(p.TDn,1);dnsv.TD = zeros(p.TDn,1);
nprior.TD = zeros(p.TDn,1);dnprior.TD = zeros(p.TDn,1);
for idx = 1:p.TDn
    [pos.TD{idx},clock_bias.TD(idx),res.TD{idx},GDOP.TD(idx),...
    nsv.TD(idx),dnsv.TD(idx),nprior.TD(idx),dnprior.TD(idx)] =...
    TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y,p.TDLambda(idx));
end
end

% [pos.MShb,clock_bias.MShb,res.MShb,GDOP.MShb] = MSHuberlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% [pos.MStk,clock_bias.MStk,res.MStk,GDOP.MStk] = MSTukeylinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% [pos.LTS,clock_bias.LTS,res.LTS,GDOP.LTS] = LTSlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% p.lambda = 4; % threshold value for TD
% [pos.TD1,clock_bias.TD1,res.TD1,GDOP.TD1] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% p.lambda = 5; % threshold value for TD
% [pos.TD2,clock_bias.TD2,res.TD2,GDOP.TD2] = TDlinear(p,xk,Pcov,H_offset,s_pos_ecef,y);
% p.lambda = 6; % threshold value for TD
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
