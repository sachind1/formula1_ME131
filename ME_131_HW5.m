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
plot(t_over_T, vt_over_Vo)
%% 2.2
% solve for beta, then get C_d and R_x
%beta = V_o * (rho*A_f*C_d/2/R_x)^1/2;

% lsqcurvefit
func_out = @(beta, t_over_T) (1/beta)*tan( (1-t_over_T) * atan(beta) );
x0 = [1]; % don't think value matters?, numel = 1 b/c only one variable (beta)
beta = lsqcurvefit(func_out, x0, t_over_T, vt_over_Vo)

% fminunc 

% C_d = (2*m*beta*atan(beta)) / (V_o*T*rho*A_f);
% R_x = (V_o*m*atan(beta)) / (beta*T);


