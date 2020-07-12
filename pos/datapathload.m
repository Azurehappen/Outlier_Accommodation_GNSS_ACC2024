function [eph_name,obs_name,IGS_name,data_base,USTEC,Grdpos,t] = datapathload(data_num,initpath)

switch(data_num)
    case 1
        % data 3
        eph_name = 'data/data200612/data200612';
        obs_name = 'data/data200612/data200612';
        IGS_name = 'data/data200612/SSRA00CAS.20C';
        data_base = [];
        USTEC = 'data/data200612/USTEC';
        Grdpos = [-2430696.646;-4704190.763;3544329.081];
        t = datetime(2020,6,11,19,26,21):seconds(1):datetime(2020,6,12,18,40,52);
    case 3
        % data 3
        eph_name = 'data/data2020/rover';
        obs_name = 'data/data2020/rover';
        IGS_name = 'data/data2020/IGS_clk10';
        % data_base = 'data/data2020/base_crfp';
        data_base = 'data/data2020/base';
        USTEC = 'data/data2020/USTEC';
        Grdpos = [-2430696.646;-4704190.763;3544329.081];
        t = datetime(2020,4,21,15,13,24):seconds(1):datetime(2020,4,22,13,29,45);
    case 4
        % data 4
        eph_name = 'data/data2020/rover';
        obs_name = 'data/data2020/rover';
        Grdpos = [-2493304.6796;-4655215.1032;3565497.5918];%JPL base station
        % load('data/data2020/base_crfp_obs.mat');
        % obs = obs_b;
        % Grdpos = [-2410446.7968;-4710490.3879;3550422.4952]; %CRFP base station
        IGS_name = 'data/data2020/IGS_clk10';
        USTEC = 'data/data2020/USTEC';
        data_base = [];
    case 6
        % data 6 Compare IGS cource
        eph_name = 'data/JPL2020144/cit11440';
        obs_name = 'data/JPL2020144/cit11440';
%         IGS_name = 'data/JPL2020144/CLK10_DREF911440.20C'; % CLK10 Germany
%         IGS_name = 'data/JPL2020144/SSRA00GFZ01440.20C'; % GFZ Germany
        IGS_name = 'data/JPL2020144/SSRA00WHU01440.20C'; % WHU
        data_base = [];
        USTEC = 'data/JPL2020144/USTEC';
        Grdpos = [-2491490.2616;-4660803.2317;3559129.0005];% Caltech base station
    case 7
        % Data 7 12 hour JPL
        eph_name = 'data/JPL2020147/JPLM00mix';
%         obs_name = 'data/JPL2020147/jplgal1470'; % GAL
        obs_name = 'data/JPL2020147/jplgps1470'; % GPS
        % obs_name= 'data/JPL2020147/jplglo1470'; % GLO
        % obs_name= 'data/JPL2020147/jplbds1470'; % BDS
        IGS_name = 'data/JPL2020147/SSRA00WHU01470.20C'; % WHU
        % IGS_name = 'data/JPL2020147/CLK10_DREF911470.20C';
        USTEC = 'data/JPL2020147/USTEC';
        data_base = [];
        Grdpos = [-2493304.6796;-4655215.1032;3565497.5918];%JPL base station
        t = datetime(2020,5,25,17,00,00):seconds(30):datetime(2020,5,26,16,59,30);
    case 8
        % Data 8, test BDS and GAL PPP
        eph_name = 'data/data202006/bdsgal';
        obs_name = 'data/data202006/bdsgal';
        IGS_name = 'data/data202006/WHU01580.20C';
        USTEC = 'data/data202006/USTEC';
        data_base = [];
        Grdpos = [-2430696.646;-4704190.763;3544329.081];
        t = datetime(2020,6,6,14,14,32):seconds(1):datetime(2020,6,6,20,41,34);
    case 9
        % Data 9, GPS,GAL PPP
        eph_name = 'data/JPL2020158/JPLM00';
        obs_name = 'data/JPL2020158/JPLMCcode'; % GPS GLO: C code
        % obs_name = 'data/JPL2020158/JPLMPcode.obs'; % GPS GLO: P code
        IGS_name = 'data/JPL2020158/WHU01580.20C';
        % IGS_name = 'data/JPL2020158/SSRA00CNE01580.20C';
        USTEC = 'data/JPL2020158/USTEC';
        data_base = [];
        Grdpos = [-2493304.6796;-4655215.1032;3565497.5918];%JPL base station
        t = datetime(2020,6,5,17,0,0):seconds(30):datetime(2020,6,7,16,59,30);
        
    case 10 % CalTrans test1
        eph_name = 'data/dataubxppp/ppp_test1';
        obs_name = 'data/dataubxppp/ppp_test1';
        IGS_name = 'data/dataubxppp/SSRA00WHU01880.20C';
        data_base = [];
        USTEC = 'data/dataubxppp/USTEC';
        %Grdpos = [-2430696.646;-4704190.763;3544329.081];
        Grdpos = [-2430697.590;-4704189.060;3544329.040];
        t = datetime(2020,7,6,2,9,52.990):seconds(1):datetime(2020,7,6,14,6,34.009); %UTC time
        
    case 11
        eph_name = 'data/dataubxppp/ppp_test1';
        obs_name = 'data/dataubxppp/jplm1880';
        IGS_name = 'data/dataubxppp/SSRA00WHU01880.20C';
        data_base = [];
        USTEC = 'data/dataubxppp/USTEC';
        Grdpos = [-2493304.6796;-4655215.1032;3565497.5918];%JPL base station
        t = datetime(2020,7,5,17,0,0):seconds(30):datetime(2020,7,6,16,59,30); 
    case 12 % JPL CalTrans test
        eph_name = 'data/CaltransTest/BRDC00WRD';
        obs_name = 'data/CaltransTest/jplm1850'; % GPS GLO: C code
        IGS_name = 'data/CaltransTest/SSRA00WHU01850.20C';
        USTEC = 'data/CaltransTest/USTEC';
        data_base = [];
        Grdpos = [-2493304.6796;-4655215.1032;3565497.5918];%JPL base station
        t = datetime(2020,7,2,17,0,0):seconds(30):datetime(2020,7,3,16,59,30);
       
end
eph_name = [initpath,eph_name];
obs_name = [initpath,obs_name];
IGS_name = [initpath,IGS_name];
USTEC = [initpath,USTEC];
if ~isempty(data_base)
    data_base = [initpath,data_base];
end
end