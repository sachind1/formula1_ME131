%% 2.1
load('straightLineTest.mat')
%rho = 1.225; % kg/m3
%m = 1755; % kg
g = 9.81; % m/s2
A_f = 1.6 + 0.00056*(m - 765); % m2
V_o = Vx(1);
T = t(end); % assuming given data goes until car effectively stops

vt_over_Vo = Vx./V_o;
t_over_T = t./T;
plot(t_over_T, vt_over_Vo) %Deliverable
%% 2.2
% solve for beta, then get C_d and R_x
%beta = V_o * (rho*A_f*C_d/2/R_x)^1/2;

% lsqcurvefit
func_out_lsq = @(beta, t_over_T) (1/beta)*tan( (1-t_over_T) * atan(beta) );
x0 = [1]; % don't think value matters?, numel = 1 b/c only one variable (beta)
beta = lsqcurvefit(func_out_lsq, x0, t_over_T, vt_over_Vo)
hold on
plot(t_over_T, (1/beta)*tan( (1-t_over_T) * atan(beta) ))

% % fminunc 
% func_out_fmu = @(beta) (1/beta)*tan( (1-blah) * atan(beta) );

C_d = (2*m*beta*atan(beta)) / (V_o*T*rho*A_f) %units?
R_x = (V_o*m*atan(beta)) / (beta*T)
%% 4.1
load('upNdownTest.mat');
% dt is sampling rate [0.01s], v1 is velocity uphill, v2 is velocity
% downhill
a1 = [];
for jdx = 2:numel(v1)
    a1 = [a1; (v1(jdx) - v1(jdx-1))/dt];
end

a2 = [];
for kdx = 2:numel(v2)
    a2 = [a2; (v2(kdx)-v2(kdx-1))/dt];
end

a1 = [a1; a1(end)];
a2 = [a2; a2(end)]; %pad

Y1 = a1;
X1 = [v1.^2 ones(numel(v1),1) -ones(numel(v1),1)]; 

Y2 = a2;
X2 = [v2.^2 ones(numel(v2),1) ones(numel(v2),1)]; % w3 is (+) b/c downhill


w1 = pinv(X1)*Y1;
w2 = pinv(X2)*Y2;

%%
Y = [a1; a2];
X = [v1.^2 ones(numel(v1),1) zeros(numel(v1),1); v2.^2 ones(numel(v2),1) zeros(numel(v2),1)]; % w3 is (+) b/c downhill

w = pinv(X)*Y

C_d_calc = w(1)*-2*m*rho*A_f
R_x_calc = w(2)*-m
theta_calc = asin(w(3)/g)









