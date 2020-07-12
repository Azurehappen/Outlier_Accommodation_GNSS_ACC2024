L_clk = 16889;
L_obt = 1406;
IGS.clk_iTOW = zeros(1,L_clk);
IGS.orbit_iTOW = zeros(1,L_obt);
IGS.GPS.clk_corr = zeros(32,L_clk);
IGS.GPS.clk_vel = zeros(32,L_clk);
IGS.GPS.clk_acc = zeros(32,L_clk);
IGS.GPS.clk_IDOC = zeros(32,L_clk);
IGS.GPS.orbit_x = zeros(32,L_obt);
IGS.GPS.orbit_y = zeros(32,L_obt);
IGS.GPS.orbit_z = zeros(32,L_obt);
IGS.GPS.orbit_xv = zeros(32,L_obt);
IGS.GPS.orbit_yv = zeros(32,L_obt);
IGS.GPS.orbit_zv = zeros(32,L_obt);
IGS.GPS.orbit_IDOC = zeros(32,L_obt);
for i = 1:L_clk
    IGS.clk_iTOW(i) = log.clock(i).iTOW;
    IGS.GPS.clk_corr(:,i) = log.clock(i).correction(1,:)';
    IGS.GPS.clk_vel(:,i) = log.clock(i).correction(2,:)';
    IGS.GPS.clk_acc(:,i) = log.clock(i).correction(3,:)';
    IGS.GPS.clk_IDOC(:,i) = log.clock(i).IOD';
end

for i = 1:L_obt
    IGS.orbit_iTOW(i) = log.orbit(i).iTOW;
    IGS.GPS.orbit_x(:,i) = log.orbit(i).correction(1,:)';
    IGS.GPS.orbit_y(:,i) = log.orbit(i).correction(2,:)';
    IGS.GPS.orbit_z(:,i) = log.orbit(i).correction(3,:)';
    IGS.GPS.orbit_xv(:,i) = log.orbit(i).correction(4,:)';
    IGS.GPS.orbit_yv(:,i) = log.orbit(i).correction(5,:)';
    IGS.GPS.orbit_zv(:,i) =log.orbit(i).correction(6,:)';
    IGS.GPS.orbit_IDOC(:,i) = log.orbit(i).IOD';
end