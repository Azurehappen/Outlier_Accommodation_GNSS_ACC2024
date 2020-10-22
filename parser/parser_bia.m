function bia = parser_bia(biapath)
% Parse code bias data from OSR bia file to matlab data file.
% Currently only parse Galileo.
% Other system can be obtained from CNES SSR data
biafile = fopen(biapath);
while (true)
    line = fgetl(biafile);
    if contains(line,'+BIAS/SOLUTION')
        line = fgetl(biafile);
        break;
    end
end
PRN_i = strfind(line,'PRN');
type_i = strfind(line,'OBS1');
val_i = strfind(line,'VALUE____');
while ~feof(biafile)
    line = fgetl(biafile);
    if strcmp(line(PRN_i),'E')
        prn = str2double(line(PRN_i+1:PRN_i+2));
        type = line(type_i:type_i+3);
        if contains(type,'C1C')
            bia.GAL.bia_C1C(prn,1) = str2double(line(val_i:val_i+8));
        end
        if contains(type,'C7Q')
            bia.GAL.bia_C7Q(prn,1) = str2double(line(val_i:val_i+8));
        end
    end
end