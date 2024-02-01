function grd = readUbxGt(file_path)

tb = readtable(file_path);

grd.t = tb.iTOW;
grd.pos = [tb.X';tb.Y';tb.Z'];
rtk_inds = tb.CarrierRangeStatus;

ind_del = rtk_inds~= 2;
grd.t(ind_del) = [];
grd.pos(:,ind_del) = [];
