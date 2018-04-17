function v_B_data = get_v(sig, Xresampled, Yresampled, Zresampled)
Ryaw = @(yaw) [
    cos(yaw), -sin(yaw), 0;
    sin(yaw),  cos(yaw), 0;
    0,         0, 1];

Rpitch = @(pitch) [
    cos(pitch), 0, sin(pitch);
    0, 1,          0;
    -sin(pitch), 0, cos(pitch)];

Rroll = @(roll) [
    1,         0,          0;
    0, cos(roll), -sin(roll);
    0, sin(roll),  cos(roll)];

Q_BtoII = @(yaw,pitch,roll) Ryaw(yaw)*Rpitch(pitch)*Rroll(roll);

yaw_data=(sig{3}.Data-sig{3}.Data(1));
pitch_data=(sig{2}.Data-sig{2}.Data(1));
roll_data=(sig{1}.Data-sig{1}.Data(1));

vx_I_signal=ts_derivative(Xresampled);
vy_I_signal=ts_derivative(Yresampled);
vz_I_signal=ts_derivative(Zresampled);

vx_I_data=vx_I_signal.Data;
vy_I_data=vy_I_signal.Data;
vz_I_data=vz_I_signal.Data;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Convert vehicle speed from inertial frame to body frame
for k=1:length(vx_I_data)
    Q1c = Q_BtoII(yaw_data(k),pitch_data(k),roll_data(k));
    out=Q1c'*[vx_I_data(k);vy_I_data(k);vz_I_data(k)];
    vx_B_data(k,1)=out(1);
    vy_B_data(k,1)=out(2);
    vz_B_data(k,1)=out(3);
end
tSpan = 86400*(sig{1}.Time-sig{1}.Time(1));

v_B_data = timeseries(sqrt(vx_B_data.^2 + vy_B_data.^2), tSpan);

end