function obs = parser_VNDLOG(path)
serverlog = fopen(path);
No_obs = 10e4;
obs = obsstctinit(No_obs);
count = 0;
while ~feof(serverlog)
    line = fgetl(serverlog);
    if contains(line,'current GPS')
        count = count + 1;
        obs.tr_prime(:,count) = sscanf(line(19:end),'%d');
        [obs.tr_week(1,count),~,obs.tr_sow(1,count)] = date2gnsst(obs.tr_prime(:,count)');
        while ~feof(serverlog)
            line = fgetl(serverlog);
            if (contains(line,'Eph_diff:')&&(line(1)=='G'))
                val = strsplit(line(2:end));
                prn = str2double(val{1});
                obs.GPS(1).data.P(prn,count) = str2double(val{6});
                obs.GPS(2).data.P(prn,count) = str2double(val{9});
            end
            if (contains(line,'Eph_diff:')&&(line(1)=='E'))
                val = strsplit(line(2:end));
                prn = str2double(val{1});
                obs.GAL(1).data.P(prn,count) = str2double(val{6});
                obs.GAL(2).data.P(prn,count) = str2double(val{9});
            end
            if (contains(line,'Eph_diff:')&&(line(1)=='C'))
                val = strsplit(line(2:end));
                prn = str2double(val{1});
                obs.BDS(1).data.P(prn,count) = str2double(val{6});
                obs.BDS(2).data.P(prn,count) = str2double(val{9});
            end
            if contains(line,'Running idx')
                break;
            end
        end
    end
end
obs = simplyobs(obs,count);
fclose(serverlog);
end