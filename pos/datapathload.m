function [eph_name,obs_name,IGS_name,data_base,code_bia,Grdpos,USTEC_folderpath] = datapathload(data_num,initpath)

switch(data_num)
    case 1
        % data 3
        eph_name = 'PPP/Multi-GNSS_test.nav';
        obs_name = 'PPP/Multi-GNSS_test';
        IGS_name = 'PPP/SSRA00CNE02940.20C';
        USTEC_folderpath = 'PPP/USTEC/';
        code_bia = 'PPP/CAS0MGXRAP_2020292_OSB.BIA'; % MGEX code bias product for GAL I channel
        data_base = [];
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;
    case 3
        eph_name = 'DGNSS_Moving/COM3_200712_013546_run1.nav';
        obs_name = 'DGNSS_Moving/COM3_200712_013546_run1';
        IGS_name = [];
        USTEC_folderpath = [];
        code_bia = [];
        data_base = 'DGNSS_Moving/rbst0711base';
        M  = readmatrix([initpath,'DGNSS_Moving/groundtruth.csv']);
        Grdpos.pos = M(:,3:5)';
        Grdpos.t = M(:,2);
    case 5
        % test CNE data PPP
        eph_name = 'VND_CNE2/COM8_0222_VND2.nav';
        obs_name = 'VND_CNE2/COM8_0222_VND2';
        IGS_name = 'VND_CNE2/SSRA00CNE00530.21C';
        USTEC_folderpath = 'VND_CNE2/USTEC/';
        code_bia = 'VND_CNE2/CAS0MGXRAP_20210490000.BIA';
        data_base = [];
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;
    case 6
        % test CNE data PPP
        eph_name = 'VND_CNE/JPLM00USA_051.rnx';
        obs_name = 'VND_CNE/COM8_0221_VND';
        IGS_name = 'VND_CNE/SSRA00CNE00520.21C';
        USTEC_folderpath = 'VND_CNE/USTEC/';
        code_bia = 'VND_CNE/CAS0MGXRAP_20210490000.BIA';
        data_base = 'VND_CNE/vndlog';
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;
    case 7
        % test CNE data PPP
        eph_name = 'VND_0227/COM8_210227_VND.nav';
        obs_name = 'VND_0227/COM8_210227_VND';
        IGS_name = [];
        USTEC_folderpath = [];
        code_bia = 'VND_0227/CAS0MGXRAP_20210490000.BIA';
        data_base = 'VND_0227/vndlog';
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;    
end
if ~isempty(eph_name)
    eph_name = [initpath,eph_name];
end
if ~isempty(obs_name)
obs_name = [initpath,obs_name];
end
if ~isempty(IGS_name)
IGS_name = [initpath,IGS_name];
end
if ~isempty(code_bia)
code_bia = [initpath,code_bia];
end
if ~isempty(USTEC_folderpath)
USTEC_folderpath = [initpath,USTEC_folderpath];
end
if ~isempty(data_base)
    data_base = [initpath,data_base];
end
end 