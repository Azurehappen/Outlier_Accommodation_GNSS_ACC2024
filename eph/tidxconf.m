    function [tidx,obt_idx,clk_idx,pos_tage,IGSdata,icb] = tidxconf(p,t_sv,prn,tidx,eph,sys_type)
% Confirm the possible time index in eph
% find the index in IGS correction data
% output the IGS data and ICB

pos_tage = 1;
IGSdata = [];
icb = [];
if p.post_mode == 1 && ~isempty(p.IGS) && p.IGS_enable == 1
    obt_idx = find(limit_tgps(t_sv - p.IGS.orbit_iTOW)<=90 & limit_tgps(t_sv - p.IGS.orbit_iTOW)>=0);
%     obt_idx = find(t_sv - p.IGS.orbit_iTOW>=0);
    if ~isempty(obt_idx)
        obt_idx = obt_idx(end);
    end
    clk_idx = find(limit_tgps(t_sv - p.IGS.clk_iTOW)<=90 & limit_tgps(t_sv - p.IGS.clk_iTOW)>=0);
%     clk_idx = find(t_sv - p.IGS.clk_iTOW>=0);
    if ~isempty(clk_idx)
        clk_idx = clk_idx(end);
    end
    if ~isempty(obt_idx) && ~isempty(clk_idx)
        switch sys_type
            case 'GPS'
                if length(tidx)>1
                    %ci = find(eph.IODC(prn,tidx) == p.IGS.GPS.clk_IDOC(prn,clk_idx));
                    ei = find(eph.IODE(prn,tidx) == p.IGS.GPS.orbit_IDOE(prn,obt_idx));
                    if ~isempty(ei)%&&isequal(ci,ei)
                        tidx = tidx(ei(end));
                    else
                        pos_tage = 0;
                    end
                else
                    if ~(...%p.IGS.GPS.clk_IDOC(prn,clk_idx) == eph.IODC(prn,tidx) &&...
                    p.IGS.GPS.orbit_IDOE(prn,obt_idx) == eph.IODE(prn,tidx))
                        pos_tage = 0;
                    end
                end
            case 'GAL'
                if prn>length(p.IGS.GAL.orbit_IDOE(:,1))
                    pos_tage = 0;
                else
                if length(tidx)>1
                    ei = find(eph.IODE(prn,tidx) == p.IGS.GAL.orbit_IDOE(prn,obt_idx));
                    if ~isempty(ei)
                        tidx = tidx(ei(end));
                    else
                        pos_tage = 0;
                    end
                else
                    if ~(p.IGS.GAL.orbit_IDOE(prn,obt_idx) == eph.IODE(prn,tidx))
                        pos_tage = 0;
                    end
                end
                end
            case 'GLO'
%                 dtr = limit_tgps(t_sv-eph.t_oc(prn,tidx));
%                 tidx = tidx(abs(dtr) == min(abs(dtr)));
                tidx=tidx(end);
            case 'BDS'
                if length(tidx)>1
                    ei = find(eph.IODE(prn,tidx) == p.IGS.BDS.orbit_IDOE(prn,obt_idx));
                    if ~isempty(ei)
                        tidx = tidx(ei(end));
                    else
                        pos_tage = 0;
                    end
                else
                    if ~(p.IGS.BDS.orbit_IDOE(prn,obt_idx) == eph.IODE(prn,tidx))
                        pos_tage = 0;
                    end
                end
        end
    else
        pos_tage = 0;
        
    end
%     cb_i = find(limit_tgps(t_sv - p.IGS.cdb_iTOW)>=0);
    if pos_tage == 1
        switch sys_type
            case 'GPS'
                IGSdata = p.IGS.GPS;
                icb = p.code_bia.GPS.bia_C1C;
%                 if ~isempty(cb_i)
%                     icb = p.IGS.GPS.code_bias_C1C(:,cb_i(end));
%                 else
%                     icb = zeros(50,1);
%                 end
            case 'GLO'
                IGSdata = p.IGS.GLO;
                icb = p.code_bia.GLO.bia_C1C;
            case 'GAL'
                IGSdata = p.IGS.GAL;
                icb = p.code_bia.GAL.bia_C1C;
            case 'BDS'
                IGSdata = p.IGS.BDS;
                icb = p.code_bia.BDS.bia_C2I;
%                 if ~isempty(cb_i)
%                     icb = p.IGS.BDS.code_bias_C2I(:,cb_i(end));
%                 else
%                     icb = zeros(50,1);
%                 end
        end
        if isempty(IGSdata)
            pos_tage = 0;
        end
    end
else
    %dtr = limit_tgps(t_sv-eph.t_oc{prn}(tidx));
    %tidx = tidx(abs(dtr) == min(abs(dtr)));
    tidx=tidx(end);
    obt_idx = [];
    clk_idx = [];
    pos_tage = 1;
end

end