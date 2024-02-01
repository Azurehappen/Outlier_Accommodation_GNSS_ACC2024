function files = dataPathLoader(data_num)

switch(data_num)
     case 1
        files.eph = 'data/ucr_ppp_moving/BRDC00WRD_S_2023237.rnx';
        files.obs = 'data/ucr_ppp_moving/rtk2.obs';
        files.ssr = 'data/ucr_ppp_moving/SSRA00WHU02370.23C';
        files.vtec = 'data/ucr_ppp_moving/SSRA00CNE02370.23C';
        files.code_bias = 'data/ucr_ppp_moving/CAS0MGXRAP_2023238_OSB.BIA';
        files.ustec_data = 'data/ucr_ppp_moving/ustec_data/';
        files.Grdpos = readUbxGt('data/ucr_ppp_moving//groud_truth_rtk.csv');
        files.preload = 'data/ucr_ppp_moving/preload.mat';
end

end 