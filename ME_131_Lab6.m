%% init
m = 5; % kg
g = 9.8; % m/s2
rho = 1.225; % kg/m3
A_f = 0.15 * 0.25 * 0.9; % m2, provided by NSCEP empirical formula
V_init = 0.5; %?
C_d = 9.0138;
R_x = 8.2029;

%%
A = 1/m;
B = rho*C_d*A_f/m;
% C = 
% D = 
Q = 1E4;
R = 1;
K = lqr(A,B,Q,R);





%%
L = 50;
wn2 = 9.8/L;
A = [0,1; wn2,0];
B = [0;-1];
C = [1,0];
Q = C'*C;
R = 1;
K1 = lqr(A,B,Q,R);
dyn = @(t,x) (A-B*K1)*x;
[tODE, xODE] = ode45(dyn, [0 15], [pi/3;0]);
subplot(3,2,1)
plot(tODE, xODE)
legend('\theta', '\theta_{dot}')
title('Q = 1, R=1')
subplot(3,2,2)
plot(tODE, K1*xODE')
ylabel('u')
Q = Q*0.01;
R = 1;
K2 = lqr(A,B,Q,R);
dyn = @(t,x) (A-B*K2)*x;
[tODE, xODE] = ode45(dyn, [0 15], [pi/3;0]);
subplot(3,2,3)
plot(tODE, xODE)
legend('\theta', '\theta_{dot}')
title('Q = 0.01, R=1')
subplot(3,2,4)
plot(tODE, K2*xODE')
ylabel('u')
Q = 10*C'*C;
R = 0.001;
K3 = lqr(A,B,Q,R);
dyn = @(t,x) (A-B*K3)*x;
[tODE, xODE] = ode45(dyn, [0 15], [pi/3;0]);
subplot(3,2,5)
plot(tODE, xODE)
legend('\theta', '\theta_{dot}')
title('Q = 10, R=0.001')
subplot(3,2,6)
plot(tODE, K3*xODE')
ylabel('u')
