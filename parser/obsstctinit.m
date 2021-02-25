function obs = obsstctinit(num)
MAXPRNGPS = 32; MAXPRNGLO = 27;
MAXPRNGAL = 36; MAXPRNBDS = 64;
% Initialize obs struct
    obs.GPS = struct;obs.GAL = struct;obs.GLO = struct;obs.BDS = struct;
    obs.GPS(1).type = 'L1 C/A';   obs.GPS(2).type = 'L2';
    obs.GPS(1).data.P = zeros(MAXPRNGPS,num);obs.GPS(2).data.P = zeros(MAXPRNGPS,num);
    obs.GPS(1).data.C = zeros(MAXPRNGPS,num);obs.GPS(2).data.C = zeros(MAXPRNGPS,num);
    obs.GPS(1).data.D = zeros(MAXPRNGPS,num);obs.GPS(2).data.D = zeros(MAXPRNGPS,num);
    obs.GPS(1).data.S = zeros(MAXPRNGPS,num);obs.GPS(2).data.S = zeros(MAXPRNGPS,num);
    obs.GPS(3).type = 'L1 P';     obs.GPS(4).type = 'L2 P';
    obs.GPS(3).data.P = zeros(MAXPRNGPS,num);obs.GPS(4).data.P = zeros(MAXPRNGPS,num);
    obs.GPS(3).data.C = zeros(MAXPRNGPS,num);obs.GPS(4).data.C = zeros(MAXPRNGPS,num);
    obs.GPS(3).data.D = zeros(MAXPRNGPS,num);obs.GPS(4).data.D = zeros(MAXPRNGPS,num);
    obs.GPS(3).data.S = zeros(MAXPRNGPS,num);obs.GPS(4).data.S = zeros(MAXPRNGPS,num);
    obs.GLO(1).type = 'L1 C/A';   obs.GLO(2).type = 'L2 C/A';
    obs.GLO(1).data.P = zeros(MAXPRNGLO,num);obs.GLO(2).data.P = zeros(MAXPRNGLO,num);
    obs.GLO(1).data.C = zeros(MAXPRNGLO,num);obs.GLO(2).data.C = zeros(MAXPRNGLO,num);
    obs.GLO(1).data.D = zeros(MAXPRNGLO,num);obs.GLO(2).data.D = zeros(MAXPRNGLO,num);
    obs.GLO(1).data.S = zeros(MAXPRNGLO,num);obs.GLO(2).data.S = zeros(MAXPRNGLO,num);
    obs.GLO(3).type = 'L1 P';     obs.GLO(4).type = 'L2 P';
    obs.GLO(3).data.P = zeros(MAXPRNGLO,num);obs.GLO(4).data.P = zeros(MAXPRNGLO,num);
    obs.GLO(3).data.C = zeros(MAXPRNGLO,num);obs.GLO(4).data.C = zeros(MAXPRNGLO,num);
    obs.GLO(3).data.D = zeros(MAXPRNGLO,num);obs.GLO(4).data.D = zeros(MAXPRNGLO,num);
    obs.GLO(3).data.S = zeros(MAXPRNGLO,num);obs.GLO(4).data.S = zeros(MAXPRNGLO,num);
    obs.GAL(1).type = 'E1';       obs.GAL(2).type = 'E5b';
    obs.GAL(1).data.P = zeros(MAXPRNGAL,num);obs.GAL(2).data.P = zeros(MAXPRNGAL,num);
    obs.GAL(1).data.C = zeros(MAXPRNGAL,num);obs.GAL(2).data.C = zeros(MAXPRNGAL,num);
    obs.GAL(1).data.D = zeros(MAXPRNGAL,num);obs.GAL(2).data.D = zeros(MAXPRNGAL,num);
    obs.GAL(1).data.S = zeros(MAXPRNGAL,num);obs.GAL(2).data.S = zeros(MAXPRNGAL,num);
    obs.BDS(1).type = 'B1';       obs.BDS(2).type = 'B2b';
	obs.BDS(1).data.P = zeros(MAXPRNBDS,num);obs.BDS(2).data.P = zeros(MAXPRNBDS,num);
    obs.BDS(1).data.C = zeros(MAXPRNBDS,num);obs.BDS(2).data.C = zeros(MAXPRNBDS,num);
    obs.BDS(1).data.D = zeros(MAXPRNBDS,num);obs.BDS(2).data.D = zeros(MAXPRNBDS,num);
    obs.BDS(1).data.S = zeros(MAXPRNBDS,num);obs.BDS(2).data.S = zeros(MAXPRNBDS,num);
    
    obs.tr_prime = NaN(6,num); obs.tr_sow = NaN(1,num); obs.tr_week = NaN(1,num);
end