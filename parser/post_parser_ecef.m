[fileName, filePath] = uigetfile({'*.pos', 'Postion Files (*.pos)'}, 'Select an Postion File');
fullPath = strcat (filePath, fileName);
[~,fileName,ext] = fileparts(fullPath);
if ext ~= ".pos"
    fprintf ('File format not supported. Please input a RINEX postion (.pos) file\n');
else
    fprintf ('Loading postions...\n \n');
end

posID = fopen(fullPath);
week = [];%GPS week
tow =  [];%GPS seconds
pos_ecef = [];%Postion in ECEF frame
num_sv = [];%The number of satellite be used in postioning
i = 1;
while ~feof(posID)
    tline = fgetl(posID);
    if strcmp(tline(1),'2')
        %Replace '/' by ' ' with the purpose of using str2num to read the time and position one time
        strtindex= strfind(tline,'/');
        tline(strtindex)=' ';
        strtindex= strfind(tline,':');
        tline(strtindex)=' '; 
        data = str2num(tline(1:77));
        [gps_week, ~, tr_gps_sow] = date2gnsst(data(1:6)); %Compute GPS second
        week = [week gps_week];
        tow = [tow tr_gps_sow];
        pos_ecef = [pos_ecef data(7:9)'];
        num_sv = [num_sv data(11)];
    end 
end
pos.week = week;
pos.tow = tow;
pos.pos_ecef = pos_ecef;
pos.num_sv = num_sv;
%Save to MAT file
save ([fileName,'_pos.mat'], 'pos');
fprintf ('Position extracted...\n \n');
%clear all; 