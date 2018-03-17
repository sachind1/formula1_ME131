%% problem 2
% 2.1
load('straightLineTest.mat')
%rho = 1.225; % kg/m3
%m = 1755; % kg
g = 9.81; % m/s2
A_f = 1.6 + 0.00056*(m - 765); % m2

C_d = 
R_x = 
beta = Vx(1) * (rho*A_f*C_d/2/R_x)^1/2; % V initial = Vx(1)

vt_over_Vo = [];
T = t(end); % assuming given data goes until car effectively stops
for idx = 1:numel(t)
    vt_over_Vo = [vt_over_Vo (1/beta)*tan( (1-t(idx)/T) * atan(beta) )];
end
t_over_T = t./T;
plot(t_over_T, vt_over_Vo)

