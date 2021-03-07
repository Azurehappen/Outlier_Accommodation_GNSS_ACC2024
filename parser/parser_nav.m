function eph = parser_nav(p,navpath)
% Parse ephemeris data from RINEX file to .mat data file. 
% Supported by RINEX version 3.03
%
%%%%%-----Reference
% ftp://igs.org/pub/data/format/rinex303.pdf
% http://acc.igs.org/misc/rinex304.pdf
%%%%%-----Input
% rinex navigation file path
%
%%%%%-----Output
% Class of constellation ephemeris
%
% Author: Azurehappen
[~,~,ext] = fileparts(navpath);
if ext == ".rnx" || ext == ".obs" || ext == ".nav"
    flag = 1;
    fprintf ('Loading ephemeris...\n \n');
else
    fprintf ('File format not supported. Please input a RINEX file\n');
end
navfile = fopen(navpath);
%----------------------------------------------------%
% Initialize variables
eph.ionoParameters=[];
eph.LeapSeconds=[];
GPS.prn_avb = []; GAL.prn_avb = []; BDS.prn_avb = []; GLO.prn_avb = [];
gcount = [];ecount = [];ccount = [];rcount = [];
g = 1; e = 1; c = 1; r = 1; %Eph count
MAXGPSPRN = p.gps.num_prn; MAXGLOPRN = p.glo.num_prn;
MAXGALPRN = p.gal.num_prn; MAXBDSPRN = p.bds.num_prn;
% Create cell for t_oc, restore different size.
GPS.t_oc = cell(MAXGPSPRN,1);GLO.t_oc = cell(MAXGLOPRN,1);
GAL.t_oc = cell(MAXGALPRN,1);BDS.t_oc = cell(MAXBDSPRN,1);
% read header
fprintf ('Reading header...\n');
ionoParameters  = []; 
while (true)    
    line = fgetl(navfile);
    lineSplit = strsplit(line);
    if contains(line,'IONOSPHERIC CORR')
        if strcmp(lineSplit(1), 'GPSA')
            ionoParameters.ionoAlpha = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GPSB')
            ionoParameters.ionoBeta = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GAL')
            ionoParameters.ionoGAL = str2double(lineSplit(2:5));
        end
            eph.ionoParameters =  ionoParameters;
    elseif contains (line,'LEAP SECONDS')
        eph.LeapSeconds = str2double(lineSplit(2));
    elseif contains(line,'END OF HEADER')
        break;
    end
end
eph.ionoParameters =  ionoParameters;
fprintf ('Finished reading the header\n \n');
%----------------------------------------------------%
% read body
fprintf ('Parsing navigation message');
while ~feof(navfile)
    line = fgetl(navfile);
    %linesp = [line(1),' ',line(2:end)]; % For example, avoid 'G15' and 'G 5' has different kind of split.
    %lineSplit = strsplit(linesp);
    sys_type = line(1);
    switch sys_type
        case{'G'} % Parse GPS ephemeris data
            prn = str2double(line(2:3));
            gcount(1,g)=prn;
            GPS.prn_avb(prn,1)=1; % prn_avb=1 means this satellite available in this dataset, to aviod matrix index exceed.
            indx = sum(gcount==prn);
            tget = str2double(strsplit(line(5:23)));
            [~,~,GPS.t_oc{prn}(indx)] = date2gpst(tget); % Time of broadcast (seconds)
            data = sscanf(line(24:end),'%f');
            GPS.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            GPS.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            GPS.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');  
            GPS.IODE(prn,indx) = data(1); % Issue of Data, Ephemeris 
            GPS.C_rs(prn,indx) = data(2); % Crs (meters) 
            GPS.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            GPS.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            GPS.C_uc(prn,indx) =  data(1); % Cuc (radians)
            GPS.e(prn,indx) = data(2); % e Eccentricity 
            GPS.C_us(prn,indx) = data(3); % Cus (radians) 
            GPS.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            GPS.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of GPS week)
            GPS.C_ic(prn,indx) = data(2); % Cic (radians)
            GPS.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            GPS.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f');   
            GPS.i_0(prn,indx) =  data(1); % i0 (radians)
            GPS.C_rc(prn,indx) = data(2); % Crc (meters) 
            GPS.Omega(prn,indx) = data(3); % omega (radians)
            GPS.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f');     
            GPS.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            GPS.CodesOnL2(prn,indx) = data(2); % Codes on L2 channel 
            GPS.week_num(prn,indx) = data(3); % GPS Week # 
            GPS.L2Pflag(prn,indx) = data(4); % L2 P data flag
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            GPS.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See GPS ICD 200H Section 20.3.3.3.1.3 
            GPS.SV_health(prn,indx) = data(2); % SV health
            GPS.TGD(prn,indx) = data(3); % TGD (seconds) 
            GPS.IODC(prn,indx) = data(4); % IODC Issue of Data, Clock 
                
            data = sscanf(fgetl(navfile),'%f'); 
            GPS.trans_time(prn,indx) = data(1); % Transmission time of message 
            GPS.fit_interval(prn,indx) = data(2); % Fit Interval in hours
            g = g +1;
            
        case{'E'} % Parse GAL ephemeris data
            prn = str2double(line(2:3));
            GAL.prn_avb(prn,1)=1;
            ecount(1,e)=prn;
            indx = sum(ecount==prn);
            tget = str2double(strsplit(line(5:23)));
            [~,~,GAL.t_oc{prn}(indx)] = date2gpst(tget); % Time of broadcast (seconds)
            data = sscanf(line(24:end),'%f');
            GAL.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            GAL.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            GAL.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');   
            GAL.IODE(prn,indx) = data(1); % Issue of Data, Ephemeris 
            GAL.C_rs(prn,indx) = data(2); % Crs (meters) 
            GAL.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            GAL.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 	  
            GAL.C_uc(prn,indx) = data(1); % Cuc (radians)
            GAL.e(prn,indx) = data(2); % e Eccentricity 
            GAL.C_us(prn,indx) = data(3); % Cus (radians) 
            GAL.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            GAL.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of GAL week)
            GAL.C_ic(prn,indx) = data(2); % Cic (radians)
            GAL.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            GAL.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            GAL.i_0(prn,indx) =  data(1); % i0 (radians)
            GAL.C_rc(prn,indx) = data(2); % Crc (meters) 
            GAL.Omega(prn,indx) = data(3); % omega (radians)
            GAL.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            GAL.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            GAL.Data_source(prn,indx) = data(2); % Data sources
            GAL.week_num(prn,indx) = data(3); % GAL Week # 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            GAL.SV_acc(prn,indx) = data(1); % SISA Signal in space accuracy (meters) 
            GAL.SV_health(prn,indx) = data(2); % SV health
            GAL.BGD_E5a(prn,indx) = data(3); % BGD E5a/E1 (seconds) 
            GAL.BGD_E5b(prn,indx) = data(4); % BGD E5b/E1 (seconds)
                
            data = sscanf(fgetl(navfile),'%f'); 
            GAL.trans_time(prn,indx) = data(1); % Transmission time of message 
            e = e +1;
            
        case{'R'}
            prn = str2double(line(2:3));
            GLO.prn_avb(prn,1)=1;
            rcount(1,r)=prn;
            indx = sum(rcount==prn);
            tget = datetime(tget); % Time of broadcast, GLO time
            tget = [tget.Year,tget.Month,tget.Day,tget.Hour,tget.Minute,tget.Second];
            [~,~,GLO.t_oc{prn}(indx)] = date2gpst(tget); % Represent GLO time by GPS time
            data = sscanf(line(24:end),'%f');
            GLO.nTauN(prn,indx) = data(1); % SV clock bias (sec) (-TauN) 
            GLO.pGammaN(prn,indx) = data(2); % SV relative frequency bias (+GammaN) 
            GLO.t_of(prn,indx) = data(3); % Message frame time (tk+nd*86400) in seconds of the UTC week 
                
            data = sscanf(fgetl(navfile),'%f'); 
            GLO.X(prn,indx) = data(1)*1000; % Satellite position X (m) 
            GLO.Xdot(prn,indx) = data(2)*1000; % velocity X dot (m/sec) 
            GLO.Xacc(prn,indx) = data(3)*1000; % X acceleration (m/sec2) 
            GLO.SV_health(prn,indx) = data(4); % health (0=OK) (Bn) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            GLO.Y(prn,indx) = data(1)*1000; % Satellite position Y (m) 
            GLO.Ydot(prn,indx) = data(2)*1000; % velocity Y dot (m/sec)
            GLO.Yacc(prn,indx) = data(3)*1000; % Y acceleration (m/sec2) 
            GLO.freq(prn,indx) = data(4); % frequency number(-7...+13) (-7...+6 ICD 5.1) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            GLO.Z(prn,indx) = data(1)*1000; % Satellite position Z (m) 
            GLO.Zdot(prn,indx) = data(2)*1000; % velocity Z dot (m/sec) 
            GLO.Zacc(prn,indx) = data(3)*1000; % Z acceleration (m/sec2) 
            GLO.age(prn,indx) = data(4); % Age of operation (days)
            r = r + 1;
        case{'C'}
            prn = str2double(line(2:3));
            BDS.prn_avb(prn,1)=1;
            ccount(1,c)=prn;
            indx = sum(ccount==prn);
            tget = str2double(strsplit(line(5:23)));
            tget = datetime(tget); % Time of broadcast, BDS time
            tget = [tget.Year,tget.Month,tget.Day,tget.Hour,tget.Minute,tget.Second];
            [~,~,BDS.t_oc{prn}(indx)] = date2gpst(tget); % Represent BDS time by GPS time
            data = sscanf(line(24:end),'%f');
            BDS.a_f0(prn,indx) = data(1); % SV clock bias (seconds) 
            BDS.a_f1(prn,indx) = data(2); % SV clock drift (sec/sec) 
            BDS.a_f2(prn,indx) = data(3); % SV clock drift rate (sec/sec2) 
            
            data = sscanf(fgetl(navfile),'%f');   
            BDS.AODE(prn,indx) = data(1); % Age of Data, Ephemeris  
            BDS.C_rs(prn,indx) = data(2); % Crs (meters) 
            BDS.Delta_n(prn,indx) = data(3); % Delta n (radians/sec) 
            BDS.M_0(prn,indx) = data(4); % M0 (radians) 
                
            data = sscanf(fgetl(navfile),'%f'); 	  
            BDS.C_uc(prn,indx) = data(1); % Cuc (radians)
            BDS.e(prn,indx) = data(2); % e Eccentricity 
            BDS.C_us(prn,indx) = data(3); % Cus (radians) 
            BDS.sqrtA(prn,indx) = data(4); % sqrt(A) (sqrt(m)) 
            
            data = sscanf(fgetl(navfile),'%f'); 
            BDS.t_oe(prn,indx) = data(1); % Toe Time of Ephemeris (sec of BDS week)
            BDS.IODE(prn,indx) = mod(floor(data(1)/720),240);
            BDS.C_ic(prn,indx) = data(2); % Cic (radians)
            BDS.Omega_0(prn,indx) = data(3); % OMEGA0 (radians)
            BDS.C_is(prn,indx) = data(4); % Cis (radians)
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            BDS.i_0(prn,indx) =  data(1); % i0 (radians)
            BDS.C_rc(prn,indx) = data(2); % Crc (meters) 
            BDS.Omega(prn,indx) = data(3); % omega (radians)
            BDS.OmegaDot(prn,indx) = data(4); % OMEGA DOT (radians/sec) 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            BDS.IDOT(prn,indx) = data(1); % IDOT (radians/sec) 
            BDS.week_num(prn,indx) = data(3); % BDT Week # 
                
            data = sscanf(fgetl(navfile),'%f'); 	    
            BDS.SV_acc(prn,indx) = data(1); % SV accuracy (meters) See BDS ICD 200H Section 20.3.3.3.1.3 
            BDS.SV_health(prn,indx) = data(2); % SV health
            BDS.TGD1(prn,indx) = data(3); % TGD1   B1/B3          (seconds) 
            BDS.TGD2(prn,indx) = data(4); % TGD2   B2/B3          (seconds) 
                
            data = sscanf(fgetl(navfile),'%f'); 
            BDS.trans_time(prn,indx) = data(1); % Transmission time of message 
            BDS.ADOC(prn,indx) = data(2); % Age of Data Clock 
            c = c +1;       
    end
end
fclose(navfile);
if length(GPS.prn_avb)<MAXGPSPRN
    GPS.prn_avb = [GPS.prn_avb; zeros(MAXGPSPRN-length(GPS.prn_avb),1)];
end
if length(GAL.prn_avb)<MAXGALPRN
    GAL.prn_avb = [GAL.prn_avb; zeros(MAXGALPRN-length(GAL.prn_avb),1)];
end
if length(GLO.prn_avb)<MAXGLOPRN
    GLO.prn_avb = [GLO.prn_avb; zeros(MAXGLOPRN-length(GLO.prn_avb),1)];
end
if length(BDS.prn_avb)<MAXBDSPRN
    BDS.prn_avb = [BDS.prn_avb; zeros(MAXBDSPRN-length(BDS.prn_avb),1)];
end
%----------------------------------------------------%
eph.GPS = GPS;
eph.GAL = GAL;
eph.GLO = GLO;
eph.BDS = BDS;
%%Save to MAT file
% save ([fileName,'_nav.mat'], 'eph');

%----------------------------------------------------%
fprintf ('\n \nEphemeris loaded:\n \n');
fprintf ('Number of GPS messages parsed: %d\n', g-1);
fprintf ('Number of GLONASS messages parsed: %d\n',r-1);
fprintf ('Number of GALileo messages parsed: %d\n',e-1);
fprintf ('Number of BDS messages parsed: %d\n',c-1);

end