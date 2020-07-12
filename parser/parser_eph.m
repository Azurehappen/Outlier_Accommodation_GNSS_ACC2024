function eph = parser_eph(p,navpath)
% Parse ephemeris data from .nav file to matlab data file. 
% Supported by RINEX version 3.03
%
%%%%%-----Reference
% ftp://igs.org/pub/data/format/rinex303.pdf
%
%%%%%-----Input
% .nav file path
%
%%%%%-----Output
% Class of constellation ephemeris
%
% Author: Wang Hu & Jean-Bernard Uwineza
% [fileName, filePath] = uigetfile({'*.nav', 'Navigation Files (*.nav)'}, 'Select a Navigation File');
% navpath = strcat (filePath, fileName);
[~,~,ext] = fileparts(navpath);

if strcmp(ext,'.nav')
fprintf ('Loading ephemeris...\n \n');
EndOfHeader = 0;
navfile = fopen(navpath);
%----------------------------------------------------%
% Initialize variables
eph.DOcreate=[];
eph.ionosphericParameters=[];
eph.LeapSeconds=[];
gps.svid_avb = []; gal.svid_avb = []; bds.svid_avb = []; glo.svid_avb = [];
gcount = [];ecount = [];ccount = [];rcount = [];
g = 1; e = 1; c = 1; r = 1; %Eph count
gps.t_oc = cell(p.gps.num_prn,1);glo.t_oc = cell(p.glo.num_prn,1);
gal.t_oc = cell(p.gal.num_prn,1);bds.t_oc = cell(p.bds.num_prn,1);
% read header
fprintf ('Reading header...\n');
while (~EndOfHeader)
    ionosphericParameters  = [];
    
    line = fgetl(navfile);
    lineSplit = strsplit(line);
    
    if contains(line,'RINEX VERSION')
        Version = lineSplit(2);
        if ~strcmp(Version,'3.03')
            error('Not the correct version, should be 3')
        end  
    elseif contains(line,'DATE')
        date = lineSplit(4);
        year = str2double(date{1,1}(1:4));
        month = str2double(date{1,1}(5:6));
        day = str2double(date{1,1}(7:8));
        eph.DOcreate = [year,month,day]; % Date of file creation 
    elseif contains(line,'IONOSPHERIC CORR')
        if strcmp(lineSplit(1), 'GPSA')
            ionosphericParameters.ionoAlpha = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GPSB')
            ionosphericParameters.ionoBeta = str2double(lineSplit(2:5));
        elseif strcmp(lineSplit(1), 'GAL')
            ionosphericParameters.ionoGAL = str2double(lineSplit(2:5));
        end
            eph.ionosphericParameters =  ionosphericParameters;
    elseif contains (line,'LEAP SECONDS')
        eph.LeapSeconds = str2double(lineSplit(2));
    elseif contains(line,'END OF HEADER')
       EndOfHeader = 1;
    end
end
fprintf ('Finished reading the header\n \n');
%----------------------------------------------------%
% read body
fprintf ('Parsing navigation message');
while ~feof(navfile)
    line = fgetl(navfile);
    linesp = [line(1),' ',line(2:end)]; % For example, avoid 'G15' and 'G 5' has different kind of split.
    lineSplit = strsplit(linesp);
    sat_id = line(1);
    switch sat_id
        case{'G'} % Parse GPS ephemeris data
            svid = str2double(line(2:3));
            gcount(1,g)=svid;
            gps.svid_avb(svid,1)=1; %gps.svid=1 means this satellite available in this dataset, to aviod matrix index exceed.
            indx = sum(gcount==svid);
            [~,~,gps.t_oc{svid}(indx)] = date2gnsst(str2double(lineSplit(3:8))); % Time of broadcast (seconds)
            gps.a_f0(svid,indx) = str2double(lineSplit(end-2)); % SV clock bias (seconds) 
            gps.a_f1(svid,indx) = str2double(lineSplit(end-1)); % SV clock drift (sec/sec) 
            gps.a_f2(svid,indx) = str2double(lineSplit(end)); % SV clock drift rate (sec/sec2) 
            
            lineSplit = strsplit(fgetl(navfile));  
            gps.IODE(svid,indx) = str2double(lineSplit(2)); % Issue of Data, Ephemeris 
            gps.C_rs(svid,indx) = str2double(lineSplit(3)); % Crs (meters) 
            gps.Delta_n(svid,indx) = str2double(lineSplit(4)); % Delta n (radians/sec) 
            gps.M_0(svid,indx) = str2double(lineSplit(5)); % M0 (radians) 
                
            lineSplit = strsplit(fgetl(navfile));	  
            gps.C_uc(svid,indx) = str2double(lineSplit(2)); % Cuc (radians)
            gps.e(svid,indx) = str2double(lineSplit(3)); % e Eccentricity 
            gps.C_us(svid,indx) = str2double(lineSplit(4)); % Cus (radians) 
            gps.sqrtA(svid,indx) = str2double(lineSplit(5)); % sqrt(A) (sqrt(m)) 
            
            lineSplit = strsplit(fgetl(navfile));
            gps.t_oe(svid,indx) = str2double(lineSplit(2)); % Toe Time of Ephemeris (sec of GPS week)
            gps.C_ic(svid,indx) = str2double(lineSplit(3)); % Cic (radians)
            gps.Omega_0(svid,indx) = str2double(lineSplit(4)); % OMEGA0 (radians)
            gps.C_is(svid,indx) = str2double(lineSplit(5)); % Cis (radians)
                
            lineSplit = strsplit(fgetl(navfile));	    
            gps.i_0(svid,indx) =  str2double(lineSplit(2)); % i0 (radians)
            gps.C_rc(svid,indx) = str2double(lineSplit(3)); % Crc (meters) 
            gps.Omega(svid,indx) = str2double(lineSplit(4)); % omega (radians)
            gps.OmegaDot(svid,indx) = str2double(lineSplit(5)); % OMEGA DOT (radians/sec) 
                
            lineSplit = strsplit(fgetl(navfile));	    
            gps.IDOT(svid,indx) = str2double(lineSplit(2)); % IDOT (radians/sec) 
            gps.CodesOnL2(svid,indx) = str2double(lineSplit(3)); % Codes on L2 channel 
            gps.week_num(svid,indx) = str2double(lineSplit(4)); % GPS Week # 
            gps.L2Pflag(svid,indx) = str2double(lineSplit(5)); % L2 P data flag
                
            lineSplit = strsplit(fgetl(navfile));	    
            gps.SV_acc(svid,indx) = str2double(lineSplit(2)); % SV accuracy (meters) See GPS ICD 200H Section 20.3.3.3.1.3 
            gps.SV_health(svid,indx) = str2double(lineSplit(3)); % SV health
            gps.TGD(svid,indx) = str2double(lineSplit(4)); % TGD (seconds) 
            gps.IODC(svid,indx) = str2double(lineSplit(5)); % IODC Issue of Data, Clock 
                
            lineSplit = strsplit(fgetl(navfile));
            gps.trans_time(svid,indx) = str2double(lineSplit(2)); % Transmission time of message 
            gps.fit_interval(svid,indx) = str2double(lineSplit(3)); % Fit Interval in hours
            g = g +1;
            
        case{'E'} % Parse GAL ephemeris data
            svid = str2double(line(2:3));
            gal.svid_avb(svid,1)=1;
            ecount(1,e)=svid;
            indx = sum(ecount==svid);
            [~,~,gal.t_oc{svid}(indx)] = date2gnsst(str2double(lineSplit(3:8))); % Time of broadcast (seconds)
            gal.a_f0(svid,indx) = str2double(lineSplit(end-2)); % SV clock bias (seconds) 
            gal.a_f1(svid,indx) = str2double(lineSplit(end-1)); % SV clock drift (sec/sec) 
            gal.a_f2(svid,indx) = str2double(lineSplit(end)); % SV clock drift rate (sec/sec2) 
            
            lineSplit = strsplit(fgetl(navfile));  
            gal.IODE(svid,indx) = str2double(lineSplit(2)); % Issue of Data, Ephemeris 
            gal.C_rs(svid,indx) = str2double(lineSplit(3)); % Crs (meters) 
            gal.Delta_n(svid,indx) = str2double(lineSplit(4)); % Delta n (radians/sec) 
            gal.M_0(svid,indx) = str2double(lineSplit(5)); % M0 (radians) 
                
            lineSplit = strsplit(fgetl(navfile));	  
            gal.C_uc(svid,indx) = str2double(lineSplit(2)); % Cuc (radians)
            gal.e(svid,indx) = str2double(lineSplit(3)); % e Eccentricity 
            gal.C_us(svid,indx) = str2double(lineSplit(4)); % Cus (radians) 
            gal.sqrtA(svid,indx) = str2double(lineSplit(5)); % sqrt(A) (sqrt(m)) 
            
            lineSplit = strsplit(fgetl(navfile));
            gal.t_oe(svid,indx) = str2double(lineSplit(2)); % Toe Time of Ephemeris (sec of gal week)
            gal.C_ic(svid,indx) = str2double(lineSplit(3)); % Cic (radians)
            gal.Omega_0(svid,indx) = str2double(lineSplit(4)); % OMEGA0 (radians)
            gal.C_is(svid,indx) = str2double(lineSplit(5)); % Cis (radians)
                
            lineSplit = strsplit(fgetl(navfile));	    
            gal.i_0(svid,indx) =  str2double(lineSplit(2)); % i0 (radians)
            gal.C_rc(svid,indx) = str2double(lineSplit(3)); % Crc (meters) 
            gal.Omega(svid,indx) = str2double(lineSplit(4)); % omega (radians)
            gal.OmegaDot(svid,indx) = str2double(lineSplit(5)); % OMEGA DOT (radians/sec) 
                
            lineSplit = strsplit(fgetl(navfile));	    
            gal.IDOT(svid,indx) = str2double(lineSplit(2)); % IDOT (radians/sec) 
            gal.Data_source(svid,indx) = str2double(lineSplit(3)); % Data sources
            gal.week_num(svid,indx) = str2double(lineSplit(4)); % GAL Week # 
                
            lineSplit = strsplit(fgetl(navfile));	    
            gal.SV_acc(svid,indx) = str2double(lineSplit(2)); % SISA Signal in space accuracy (meters) 
            gal.SV_health(svid,indx) = str2double(lineSplit(3)); % SV health
            gal.BGD_E5a(svid,indx) = str2double(lineSplit(4)); % BGD E5a/E1 (seconds) 
            gal.BGD_E5B(svid,indx) = str2double(lineSplit(5)); % BGD E5b/E1 (seconds)
                
            lineSplit = strsplit(fgetl(navfile));
            gal.trans_time(svid,indx) = str2double(lineSplit(2)); % Transmission time of message 
            e = e +1;
            
        case{'R'}
            svid = str2double(line(2:3));
            glo.svid_avb(svid,1)=1;
            rcount(1,r)=svid;
            indx = sum(rcount==svid);
            [~,~,glo.t_oc{svid}(indx)] = date2gnsst(str2double(lineSplit(3:8))); % Time of broadcast
            glo.nTauN(svid,indx) = str2double(lineSplit(end-2)); % SV clock bias (sec) (-TauN) 
            glo.pGammaN(svid,indx) = str2double(lineSplit(end-1)); % SV relative frequency bias (+GammaN) 
            glo.t_of(svid,indx) = str2double(lineSplit(end)); % Message frame time (tk+nd*86400) in seconds of the UTC week 
                
            lineSplit = strsplit(fgetl(navfile));
            glo.X(svid,indx) = str2double(lineSplit(2))*1000; % Satellite position X (m) 
            glo.Xdot(svid,indx) = str2double(lineSplit(3))*1000; % velocity X dot (m/sec) 
            glo.Xacc(svid,indx) = str2double(lineSplit(4))*1000; % X acceleration (m/sec2) 
            glo.SV_health(svid,indx) = str2double(lineSplit(5)); % health (0=OK) (Bn) 
                
            lineSplit = strsplit(fgetl(navfile));
            glo.Y(svid,indx) = str2double(lineSplit(2))*1000; % Satellite position Y (m) 
            glo.Ydot(svid,indx) = str2double(lineSplit(3))*1000; % velocity Y dot (m/sec)
            glo.Yacc(svid,indx) = str2double(lineSplit(4))*1000; % Y acceleration (m/sec2) 
            glo.freq(svid,indx) = str2double(lineSplit(5)); % frequency number(-7...+13) (-7...+6 ICD 5.1) 
                
            lineSplit = strsplit(fgetl(navfile));
            glo.Z(svid,indx) = str2double(lineSplit(2))*1000; % Satellite position Z (m) 
            glo.Zdot(svid,indx) = str2double(lineSplit(3))*1000; % velocity Z dot (m/sec) 
            glo.Zacc(svid,indx) = str2double(lineSplit(4))*1000; % Z acceleration (m/sec2) 
            glo.age(svid,indx) = str2double(lineSplit(5)); % Age of operation (days)
            r = r + 1;
        case{'C'}
            svid = str2double(line(2:3));
            bds.svid_avb(svid,1)=1;
            ccount(1,c)=svid;
            indx = sum(ccount==svid);
            [~,~,bds.t_oc{svid}(indx)] = date2gnsst(str2double(lineSplit(3:8))); % Time of broadcast
            bds.a_f0(svid,indx) = str2double(lineSplit(end-2)); % SV clock bias (seconds) 
            bds.a_f1(svid,indx) = str2double(lineSplit(end-1)); % SV clock drift (sec/sec) 
            bds.a_f2(svid,indx) = str2double(lineSplit(end)); % SV clock drift rate (sec/sec2) 
            
            lineSplit = strsplit(fgetl(navfile));  
            bds.IODE(svid,indx) = str2double(lineSplit(2)); % Age of Data, Ephemeris  
            bds.C_rs(svid,indx) = str2double(lineSplit(3)); % Crs (meters) 
            bds.Delta_n(svid,indx) = str2double(lineSplit(4)); % Delta n (radians/sec) 
            bds.M_0(svid,indx) = str2double(lineSplit(5)); % M0 (radians) 
                
            lineSplit = strsplit(fgetl(navfile));	  
            bds.C_uc(svid,indx) = str2double(lineSplit(2)); % Cuc (radians)
            bds.e(svid,indx) = str2double(lineSplit(3)); % e Eccentricity 
            bds.C_us(svid,indx) = str2double(lineSplit(4)); % Cus (radians) 
            bds.sqrtA(svid,indx) = str2double(lineSplit(5)); % sqrt(A) (sqrt(m)) 
            
            lineSplit = strsplit(fgetl(navfile));
            bds.t_oe(svid,indx) = str2double(lineSplit(2)); % Toe Time of Ephemeris (sec of bds week)
            bds.C_ic(svid,indx) = str2double(lineSplit(3)); % Cic (radians)
            bds.Omega_0(svid,indx) = str2double(lineSplit(4)); % OMEGA0 (radians)
            bds.C_is(svid,indx) = str2double(lineSplit(5)); % Cis (radians)
                
            lineSplit = strsplit(fgetl(navfile));	    
            bds.i_0(svid,indx) =  str2double(lineSplit(2)); % i0 (radians)
            bds.C_rc(svid,indx) = str2double(lineSplit(3)); % Crc (meters) 
            bds.Omega(svid,indx) = str2double(lineSplit(4)); % omega (radians)
            bds.OmegaDot(svid,indx) = str2double(lineSplit(5)); % OMEGA DOT (radians/sec) 
                
            lineSplit = strsplit(fgetl(navfile));	    
            bds.IDOT(svid,indx) = str2double(lineSplit(2)); % IDOT (radians/sec) 
            bds.week_num(svid,indx) = str2double(lineSplit(4)); % BDT Week # 
                
            lineSplit = strsplit(fgetl(navfile));	    
            bds.SV_acc(svid,indx) = str2double(lineSplit(2)); % SV accuracy (meters) See bds ICD 200H Section 20.3.3.3.1.3 
            bds.SV_health(svid,indx) = str2double(lineSplit(3)); % SV health
            bds.TGD1(svid,indx) = str2double(lineSplit(4)); % TGD1   B1/B3          (seconds) 
            bds.TGD2(svid,indx) = str2double(lineSplit(5)); % TGD2   B2/B3          (seconds) 
                
            lineSplit = strsplit(fgetl(navfile));
            bds.trans_time(svid,indx) = str2double(lineSplit(2)); % Transmission time of message 
            bds.IDOC(svid,indx) = str2double(lineSplit(3)); % Age of Data Clock 
            c = c +1;       
        otherwise
            fprintf ('Unkown constellation\n');
    end
end
fclose(navfile);
if length(gps.svid_avb)<p.gps.num_prn
    gps.svid_avb = [gps.svid_avb; zeros(33-length(gps.svid_avb),1)];
end
if length(gal.svid_avb)<p.gal.num_prn
    gal.svid_avb = [gal.svid_avb; zeros(38-length(gal.svid_avb),1)];
end
if length(glo.svid_avb)<p.glo.num_prn
    glo.svid_avb = [glo.svid_avb; zeros(34-length(glo.svid_avb),1)];
end
if length(bds.svid_avb)<p.bds.num_prn
    bds.svid_avb = [bds.svid_avb; zeros(38-length(bds.svid_avb),1)];
end
%----------------------------------------------------%
eph.gps = gps;
eph.gal=gal;
eph.glo=glo;
eph.bds=bds;
%%Save to MAT file
% save ([fileName,'_nav.mat'], 'eph');
else
fprintf ('File format not suppoerted. Please input a RINEX navigation (.nav) file\n');
end
%----------------------------------------------------%
fprintf ('\n \n Ephemeris loaded correctly\n \n');
fprintf ('Number of GPS messages parsed: %d\n', g-1);
fprintf ('Number of GLONASS messages parsed: %d\n',r-1);
fprintf ('Number of Galileo messages parsed: %d\n',e-1);
fprintf ('Number of BDS messages parsed: %d\n',c-1);

end