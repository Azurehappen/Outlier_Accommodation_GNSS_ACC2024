function IGS = parser_IGS(IGS_path)


IGSfile = fopen(IGS_path);
obt_i = 0;
clk_i = 0;
cdb_i = 0;
type = 0;
IGS.GPS=[];
IGS.GLO=[];
IGS.GAL=[];
IGS.BDS=[];
while ~feof(IGSfile)
   line = fgetl(IGSfile);
   if strcmp(line(1),'>') % new period
       linesplit = strsplit(line);
       switch linesplit{2}
           case 'ORBIT'
               obt_i = obt_i+1;
               [~,~,IGS.orbit_iTOW(:,obt_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               
               type = 1;
           case 'CLOCK'
               clk_i = clk_i+1;
%                IGS.clk_tr(:,clk_i) = str2double(linesplit(3:8))'; % [year;month;date;hour;minute;second]
               [~,~,IGS.clk_iTOW(:,clk_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               
               type = 2;
           case 'CODE_BIAS'
               cdb_i = cdb_i+1;
               [~,~,IGS.cdb_iTOW(:,cdb_i)] = date2gnsst(str2double(linesplit(3:8)));%  GPS seconds
               
               type = 3;
           case 'PHASE_BIAS'
               type = 4;
       end
   else
       linesp = [line(1),' ',line(2:end)];
       linesplit = strsplit(linesp);
       i = str2double(linesplit(2));
       switch type
           case 1
               switch linesplit{1}
                   case 'G'
                       IGS.GPS.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       IGS.GPS.orbit_x(i,obt_i) = str2double(linesplit(4));
                       IGS.GPS.orbit_y(i,obt_i) = str2double(linesplit(5));
                       IGS.GPS.orbit_z(i,obt_i) = str2double(linesplit(6));
                       IGS.GPS.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       IGS.GPS.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       IGS.GPS.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'R'
                       IGS.GLO.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       IGS.GLO.orbit_x(i,obt_i) = str2double(linesplit(4));
                       IGS.GLO.orbit_y(i,obt_i) = str2double(linesplit(5));
                       IGS.GLO.orbit_z(i,obt_i) = str2double(linesplit(6));
                       IGS.GLO.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       IGS.GLO.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       IGS.GLO.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'E'
                       IGS.GAL.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       IGS.GAL.orbit_x(i,obt_i) = str2double(linesplit(4));
                       IGS.GAL.orbit_y(i,obt_i) = str2double(linesplit(5));
                       IGS.GAL.orbit_z(i,obt_i) = str2double(linesplit(6));
                       IGS.GAL.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       IGS.GAL.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       IGS.GAL.orbit_zv(i,obt_i) = str2double(linesplit(9));
                   case 'C'
                       IGS.BDS.orbit_IDOE(i,obt_i) = str2double(linesplit(3));
                       IGS.BDS.orbit_x(i,obt_i) = str2double(linesplit(4));
                       IGS.BDS.orbit_y(i,obt_i) = str2double(linesplit(5));
                       IGS.BDS.orbit_z(i,obt_i) = str2double(linesplit(6));
                       IGS.BDS.orbit_xv(i,obt_i) = str2double(linesplit(7));
                       IGS.BDS.orbit_yv(i,obt_i) = str2double(linesplit(8));
                       IGS.BDS.orbit_zv(i,obt_i) = str2double(linesplit(9));
               end
           case 2 
               switch linesplit{1}
                   case 'G'
                       IGS.GPS.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       IGS.GPS.clk_corr(i,clk_i) = str2double(linesplit(4));
                       IGS.GPS.clk_vel(i,clk_i) = str2double(linesplit(5));
                       IGS.GPS.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'R'
                       IGS.GLO.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       IGS.GLO.clk_corr(i,clk_i) = str2double(linesplit(4));
                       IGS.GLO.clk_vel(i,clk_i) = str2double(linesplit(5));
                       IGS.GLO.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'E'
                       IGS.GAL.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       IGS.GAL.clk_corr(i,clk_i) = str2double(linesplit(4));
                       IGS.GAL.clk_vel(i,clk_i) = str2double(linesplit(5));
                       IGS.GAL.clk_acc(i,clk_i) = str2double(linesplit(6));
                   case 'C'
                       IGS.BDS.clk_IDOC(i,clk_i) = str2double(linesplit(3));
                       IGS.BDS.clk_corr(i,clk_i) = str2double(linesplit(4));
                       IGS.BDS.clk_vel(i,clk_i) = str2double(linesplit(5));
                       IGS.BDS.clk_acc(i,clk_i) = str2double(linesplit(6));
               end
           case 3
               switch linesplit{1}
                   case 'G'
                       IGS.GPS.code_bias_L1(i,cdb_i) = str2double(linesplit(5));
                   case 'R'
                       IGS.GLO.code_bias_L1(i,cdb_i) = str2double(linesplit(5));
                   case 'E'
                       IGS.GAL.code_bias_L1(i,cdb_i) = str2double(linesplit(5));
                   case 'C'
                       IGS.BDS.code_bias_L1(i,cdb_i) = str2double(linesplit(5));    
               end
           case 4
               
       end
              
   end
end    
fclose(IGSfile);
end