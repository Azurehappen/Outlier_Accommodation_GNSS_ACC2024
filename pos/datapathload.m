function [eph_name,obs_name,IGS_name,data_base,code_bia,Grdpos,USTEC_folderpath] = datapathload(data_num,initpath)

switch(data_num)
    case 1
        % data 3
        eph_name = 'Multi-GNSS_test';
        obs_name = 'Multi-GNSS_test';
        IGS_name = 'SSRA00CNE02940.20C';
        USTEC_folderpath = "data/USTEC/";
        code_bia = 'CAS0MGXRAP_2020291_OSB.BIA'; % MGEX code bias product for GAL I channel
        data_base = [];
        Grdpos = [-2430697.699;-4704189.201;3544329.103];
    case 2
        % Data 8, test BDS and GAL PPP
        eph_name = 'data202006/bdsgal';
        obs_name = 'data202006/bdsgal';
        IGS_name = 'data202006/WHU01580.20C';
        USTEC_folderpath = 'data202006/USTEC';
        code_bia = [];
        data_base = [];
        Grdpos = [-2430696.646;-4704190.763;3544329.081];
end
eph_name = [initpath,eph_name];
obs_name = [initpath,obs_name];
IGS_name = [initpath,IGS_name];
if ~isempty(data_base)
    data_base = [initpath,data_base];
end
end 