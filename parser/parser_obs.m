function obs = parser_obs(obspath)
% Parse observation data from .obs file to matlab data file. 
% Supported by RINEX version 3.03
%
%%%%%-----Reference
% ftp://igs.org/pub/data/format/rinex303.pdf
%
%%%%%-----Input
% .obs file path
%
%%%%%-----Output
% Class of constellation observatoin
%
% Specification
%           Output obs will include GPS, GAL, GLO and BDS
%           gps.P1: psedorange from f1 frequency. gps.P2: psedorange from f2 frequency. 
%           gps.C1: carrier phase  from f1 frequency. gps.C2: carrier phase from f2 frequency. 
%           gps.D1: doppler from f1 frequency. gps.D2: doppler from f2 frequency. 
%           gps.S1: signal strength from f1 frequency. gps.S2: signal strength from f2 frequency. 
%           gal.P1/C1/D1/S1 for f1 frequency, gal.P2/C2/D2/S2 for f5 frequency
%           bds.P1/C1/D1/S1 for f1 frequency, bds.P2/C2/D2/S2 for f2 frequency
%           glo.P1/C1/D1/S1 for f1 frequency, glo.P2/C2/D2/S2 for f2 frequency
%           ~.dual = 0 means single frequency, 1 means dual frequency
%           specific definition of f1 or f2 depends on the message from 'SYS / # / OBS TYPES' in obs
%           see RINEX 3.03 table A2
%
% Author: Wang Hu
[~,~,ext] = fileparts(obspath);

if strcmp(ext,'.obs')
fprintf ('Loading observations...\n \n');
obsfile = fopen(obspath);
%-----------------------------------%
% Initialization
obs.GPS = [];obs.GAL = [];obs.GLO = [];obs.BDS = [];
% read header
fprintf ('Reading header...\n');
while (true)  
    line = fgetl(obsfile);                                                   %get line    
    %     if contains(line,'APPROX POSITION XYZ')                                  % Receiver aprox position
    %         XYZ_station=real(str2doubleq(line(1:60)));
    %     end
    % find number of constellations available and put the info in obs.GNSS
    if contains(line,'SYS / # / OBS TYPES')
        constellation = line(1);
        num_obstype  = str2double(line(5:6));       % number of observables
        if constellation        == 'G'
            fprintf('File contains GPS observations \n')
            obs.GPS.num_obstype = num_obstype;
            obs.GPS.f1 = gpsfreq(str2double(line(9)));
            obs.GPS.P1 = [];
            if num_obstype > 4
                obs.GPS_dual = 1;
                obs.GPS.f2 = gpsfreq(str2double(line(25)));
            elseif num_obstype == 4
                obs.GPS_dual = 0;
            end
        elseif constellation    == 'R'
            fprintf('File contains GLO observations \n')
            obs.GLO.num_obstype = num_obstype;
            obs.GLO.f1 = glofreq(str2double(line(9)));
            obs.GLO.P1 = [];
            if num_obstype > 4
                obs.GLO_dual = 1;
            elseif num_obstype == 4
                obs.GLO_dual = 0;
                obs.GLO.f2 = glofreq(str2double(line(25)));
            end
        elseif constellation    == 'E'
            fprintf('File contains GAL observations \n')
            obs.GAL.num_obstype = num_obstype;
            obs.GAL.f1 = galfreq(str2double(line(9)));
            obs.GAL.P1 = [];
            if num_obstype > 4
                obs.GAL_dual = 1;
            elseif num_obstype == 4
                obs.GAL_dual = 0;
                obs.GAL.f2 = galfreq(str2double(line(25)));
            end
        elseif constellation    == 'C'
            fprintf('File contains BDS observations \n')
            obs.BDS.num_obstype = num_obstype;
            obs.BDS.f1 = bdsfreq(str2double(line(9)));
            obs.BDS.P1 = [];
            if num_obstype > 4
                obs.BDS_dual = 1;
            elseif num_obstype == 4
                obs.BDS_dual = 0;
                obs.BDS.f2 = bdsfreq(str2double(line(25)));
            end
        end        
    elseif contains(line,'END OF HEADER')
        break;                                                              % End of header loop
    end
end
%-----------------------------------%
% read observables
count = 0;
fprintf ('Parsing observables \n');
while ~feof(obsfile)
    line = fgetl(obsfile);
    if strcmp(line(1),'>') % new observables
        count = count+1;
        lineSplit = strsplit(line);
        obs.tr_prime(:,count) = str2double(lineSplit(2:7))'; % [year;month;date;hour;minute;second]
        [~,~,obs.tr_sow(1,count)] = date2gnsst(str2double(lineSplit(2:7))); %  GPS seconds
    else
        switch line(1)
            case{'G'}
                prn = str2double(line(2:3));
                obs.GPS.P1(prn,count)= str2double(line(6:17));
                obs.GPS.C1(prn,count)= str2double(line(21:33));
                obs.GPS.D1(prn,count)= str2double(line(41:49));
                obs.GPS.S1(prn,count)= str2double(line(60:65));
                if obs.GPS_dual == 1 % decide if it's dual frequency
                    obs.GPS.P2(prn,count)= str2double(line(70:81));
                    obs.GPS.C2(prn,count)= str2double(line(86:97));
                    obs.GPS.D2(prn,count)= str2double(line(105:113));
                    obs.GPS.S2(prn,count)= str2double(line(124:129));
                end
            case{'E'}
                prn = str2double(line(2:3));
                obs.GAL.P1(prn,count)= str2double(line(6:17));
                obs.GAL.C1(prn,count)= str2double(line(21:33));
                obs.GAL.D1(prn,count)= str2double(line(41:49));
                obs.GAL.S1(prn,count)= str2double(line(60:65));
                if obs.GAL_dual == 1
                    switch obs.GAL.num_obstype
                        case 8
                            obs.GAL.P2(prn,count)= str2double(line(70:81));
                            obs.GAL.C2(prn,count)= str2double(line(86:97));
                            obs.GAL.D2(prn,count)= str2double(line(105:113));
                            obs.GAL.S2(prn,count)= str2double(line(124:129));
                        case 12
                            obs.GAL.P2(prn,count)= str2double(line(134:145));
                            obs.GAL.C2(prn,count)= str2double(line(150:161));
                            obs.GAL.D2(prn,count)= str2double(line(169:177));
                            obs.GAL.S2(prn,count)= str2double(line(188:193));
                    end          
                end
            case{'R'}
                prn = str2double(line(2:3));
                obs.GLO.P1(prn,count)= str2double(line(6:17));
                obs.GLO.C1(prn,count)= str2double(line(21:33));
                obs.GLO.D1(prn,count)= str2double(line(41:49));
                obs.GLO.S1(prn,count)= str2double(line(60:65));
                if obs.GLO_dual == 1
                    obs.GLO.P2(prn,count)= str2double(line(70:81));
                    obs.GLO.C2(prn,count)= str2double(line(86:97));
                    obs.GLO.D2(prn,count)= str2double(line(105:113));
                    obs.GLO.S2(prn,count)= str2double(line(124:129));
                end
            case{'C'}
                prn = str2double(line(2:3));
                obs.BDS.P1(prn,count)= str2double(line(6:17));
                obs.BDS.C1(prn,count)= str2double(line(21:33));
                obs.BDS.D1(prn,count)= str2double(line(41:49));
                obs.BDS.S1(prn,count)= str2double(line(60:65));
                if obs.BDS_dual == 1
                    obs.BDS.P2(prn,count)= str2double(line(70:81));
                    obs.BDS.C2(prn,count)= str2double(line(86:97));
                    obs.BDS.D2(prn,count)= str2double(line(105:113));
                    obs.BDS.S2(prn,count)= str2double(line(124:129));
                end
        end
    end 
end
fclose(obsfile);
%%Save to MAT file
% save ([fileName,'_obs.mat'], 'obs');
else
fprintf ('File format not suppoerted. Please input a RINEX navigation (.obs) file\n');
end
%----------------------------------------------------%
fprintf ('\n \n Observables loaded correctly\n \n');

% Define the function that to parse the type of frequency
    function type = gpsfreq(num)
        switch num
            case 1
                type = 'L1';
            case 2
                type = 'L2';
            case 5
                type = 'L5';
            otherwise
                type = [];
        end
    end
    function type = galfreq(num)
        switch num
            case 1
                type = 'E1';
            case 5
                type = 'E5a';
            case 6
                type = 'E6';
            case 7
                type = 'E5b';
            case 8 
                type = 'E5ab';
            otherwise
                type = [];
        end
    end
    function type = glofreq(num)
        switch num
            case 1
                type = 'G1';
            case 2
                type = 'G2';
            otherwise
                type = [];
        end
    end
    function type = bdsfreq(num)
        switch num
            case 2
                type = 'B1';
            case 6
                type = 'B3';
            case 7
                type = 'B2';
            otherwise
                type = [];
        end
    end


end

