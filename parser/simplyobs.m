function obs = simplyobs(obs,count)
% Delete empty data in obs strcut

if count < length(obs.tr_sow)
    obs.tr_prime(:,count+1:end)=[]; 
    obs.tr_sow(:,count+1:end)=[];
    obs.tr_week(:,count+1:end)=[];
    obs.GPS(1).data.P(:,count+1:end)=[];  obs.GPS(2).data.P(:,count+1:end)=[];
    obs.GPS(1).data.C(:,count+1:end)=[];  obs.GPS(2).data.C(:,count+1:end)=[];
    obs.GPS(1).data.D(:,count+1:end)=[];  obs.GPS(2).data.D(:,count+1:end)=[];
    obs.GPS(1).data.S(:,count+1:end)=[];  obs.GPS(2).data.S(:,count+1:end)=[];
    obs.GPS(3).data.P(:,count+1:end)=[];  obs.GPS(4).data.P(:,count+1:end)=[];
    obs.GPS(3).data.C(:,count+1:end)=[];  obs.GPS(4).data.C(:,count+1:end)=[];
    obs.GPS(3).data.D(:,count+1:end)=[];  obs.GPS(4).data.D(:,count+1:end)=[];
    obs.GPS(3).data.S(:,count+1:end)=[];  obs.GPS(4).data.S(:,count+1:end)=[];
    obs.GLO(1).data.P(:,count+1:end)=[];  obs.GLO(2).data.P(:,count+1:end)=[];
    obs.GLO(1).data.C(:,count+1:end)=[];  obs.GLO(2).data.C(:,count+1:end)=[];
    obs.GLO(1).data.D(:,count+1:end)=[];  obs.GLO(2).data.D(:,count+1:end)=[];
    obs.GLO(1).data.S(:,count+1:end)=[];  obs.GLO(2).data.S(:,count+1:end)=[];
    obs.GLO(3).data.P(:,count+1:end)=[];  obs.GLO(4).data.P(:,count+1:end)=[];
    obs.GLO(3).data.C(:,count+1:end)=[];  obs.GLO(4).data.C(:,count+1:end)=[];
    obs.GLO(3).data.D(:,count+1:end)=[];  obs.GLO(4).data.D(:,count+1:end)=[];
    obs.GLO(3).data.S(:,count+1:end)=[];  obs.GLO(4).data.S(:,count+1:end)=[];
    obs.GAL(1).data.P(:,count+1:end)=[];  obs.GAL(2).data.P(:,count+1:end)=[];
    obs.GAL(1).data.C(:,count+1:end)=[];  obs.GAL(2).data.C(:,count+1:end)=[];
    obs.GAL(1).data.D(:,count+1:end)=[];  obs.GAL(2).data.D(:,count+1:end)=[];
    obs.GAL(1).data.S(:,count+1:end)=[];  obs.GAL(2).data.S(:,count+1:end)=[];
	obs.BDS(1).data.P(:,count+1:end)=[];  obs.BDS(2).data.P(:,count+1:end)=[];
    obs.BDS(1).data.C(:,count+1:end)=[];  obs.BDS(2).data.C(:,count+1:end)=[];
    obs.BDS(1).data.D(:,count+1:end)=[];  obs.BDS(2).data.D(:,count+1:end)=[];
    obs.BDS(1).data.S(:,count+1:end)=[];  obs.BDS(2).data.S(:,count+1:end)=[];
end

end