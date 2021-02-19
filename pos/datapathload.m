function [eph_name,obs_name,IGS_name,data_base,code_bia,Grdpos,USTEC_folderpath] = datapathload(data_num,initpath)

switch(data_num)
    case 1
        % data 3
        eph_name = 'Multi-GNSS_test';
        obs_name = 'Multi-GNSS_test';
        IGS_name = 'SSRA00CNE02940.20C';
        USTEC_folderpath = "data/USTEC/";
        code_bia = 'CAS0MGXRAP_2020292_OSB.BIA'; % MGEX code bias product for GAL I channel
        data_base = [];
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;
    case 2
        % Data 8, test BDS and GAL PPP
        eph_name = 'data202006/bdsgal';
        obs_name = 'data202006/bdsgal';
        IGS_name = 'data202006/WHU01580.20C';
        USTEC_folderpath = 'data202006/USTEC';
        code_bia = [];
        data_base = [];
        Grdpos.pos = [-2430697.699;-4704189.201;3544329.103];
        Grdpos.t = NaN;
    case 3
        eph_name = 'COM3_200712_013546_run1';
        obs_name = 'COM3_200712_013546_run1';
        IGS_name = [];
        USTEC_folderpath = [];
        code_bia = [];
        data_base = 'rbst0711base';
        M  = readmatrix([initpath,'groundtruth.csv']);
        Grdpos.pos = M(:,3:5)';
        Grdpos.t = M(:,2);
end
eph_name = [initpath,eph_name];
obs_name = [initpath,obs_name];
IGS_name = [initpath,IGS_name];
if ~isempty(data_base)
    data_base = [initpath,data_base];
end
end 