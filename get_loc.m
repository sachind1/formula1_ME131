function [X,Y, Z] = get_loc(sig)

latitude_data = sig{1,17}.Data-sig{1,17}.Data(1);
longitude_data = sig{1,16}.Data-sig{1,16}.Data(1);
altitude_data = sig{1,18}.Data-sig{1,18}.Data(1);
pre_center = [latitude_data(1),longitude_data(1)];

local_output = lla2flat([latitude_data, longitude_data, altitude_data], pre_center,sig{3}.Data(1)*180/pi-90,altitude_data(1));
X = local_output(:,2);
Y= local_output(:,1);
Z = local_output(:,3);

Xsig=timeseries(X,sig{17}.Time); 
Ysig=timeseries(Y,sig{17}.Time);
Zsig = timeseries(Z,sig{17}.Time);

X=resample(Xsig, sig{1}.Time,'linear');
Y=resample(Ysig, sig{1}.Time,'linear');
Z = resample(Zsig, sig{1}.Time, 'linear');

end