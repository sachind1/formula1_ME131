function [x_store,y_store] = simulate_system(sig, N, lf)

% actual x, y sampled at IMU time
[X_out,Y_out,Z_out] = get_loc(sig);
X = @(t) interp1((X_out.time - X_out.Time(1))*24*3600, X_out.data, t);
Y = @(t) interp1((Y_out.time - Y_out.Time(1))*24*3600, Y_out.data, t);

% actual v sampled at IMU time
V_out = get_v(sig, X_out, Y_out, Z_out);
V = @(t) interp1(V_out.time, V_out.data, t);

% actual yaw rate (inherently sampled at IMU time)
Yaw = @(t) interp1((sig{1,3}.Time- sig{1,3}.Time(1))*24*3600, smooth(unwrap(sig{1,3}.Data)-sig{1,3}.Data(1)),t);

% set intervals accrding to given N
imu_time = (sig{1,1}.Time - sig{1,1}.Time(1))*24*3600; 
interval_time = (imu_time(end)-imu_time(2))/N;
interval_length = floor((length(imu_time)-1) / N);

% motor and servo inputs (sampled at IMU time)
servo_pwm = resample(sig{1,11}, sig{1}.Time, 'linear');
servo_pwm = servo_pwm.Data;
u_s =@(t) interp1(imu_time, servo_pwm, t);

motor_pwm = resample(sig{1,10}, sig{1}.Time, 'linear');
motor_pwm = motor_pwm.Data;
u_m =@(t) interp1(imu_time, motor_pwm, t);

% Parameters
L = 0.32;
lr = L - lf;

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% Deliverable 1:
% fill in your a and b values to define your velocity and steering mappings
ab_motor = [-0.5742 0.009646]; % TO DO: fill in vector
ab_servo = [-0.50 0.01]; % TO DO: fill in vector

% write a function handle that outputs the steering angle applied to the
% car during the experiment as a function of time
steering_angle = @(t) ab_servo(1)*u_s(t) + ab_servo(2); %TO DO: define here

% write a function handle that outputs the beta angle during the experiment
% as a function of time
beta = @(t) atan( (lf*tan(0)+lr*tan(steering_angle(t))) / (lf+lr) ); % TO DO: define here

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% we will store simulated x and y here
x_store=[];
y_store=[];

% initialize our first simulation step
v0 = 0;
phi0 = Yaw(imu_time(1));
x0 = 0;
y0 = 0;
    
%% Deliverable 2: 
% complete the for-loop using the kinematic bicycle model
for i=1:N
    tSpan = [imu_time(2)+(i-1)*interval_time, imu_time(2)+(i)*interval_time]; %in seconds
    
    % write a function handle for the velocity dynamics of the car: 
    v_dyn = @(t, v, u_m) ab_motor(1)*v(t) + ab_motor(2)*u_m(t); % TO DO: define here
    velocity = ode45(@(t,v) v_dyn(t,v,u_m), tSpan, v0);
    v = @(t) deval(velocity,t);
    
    % write a function handle for the yaw angle dynamics of the car:   PSI
    phi_dyn = @(t) v(t)*cos(beta(t)) * (tan(steering_angle(t)) - tan(0)) / (lf+lr); % TO DO: define here
    phi_solved = ode45(@(t,v) phi_dyn(t), tSpan, phi0);
    phi = @(t) deval(phi_solved,t);
    
    % write a function handle for the x position dynamics of the car:
    x_dyn = @(t) v(t)*cos(phi(t) + beta(t)); % TO DO: define here
    x_solved = ode45(@(t,v) x_dyn(t), tSpan, x0);
    x = @(t) deval(x_solved, t);
    
    % write a function handle for the y position dynamics of the car:
    y_dyn = @(t) v(t)*sin(phi(t) + beta(t)); % TO DO: define here
    y_solved = ode45(@(t,v) y_dyn(t), tSpan, y0);
    y = @(t) deval(y_solved, t);
    
    % simulate x and y values forward using the kinematic bicycle model
    x_store = [x_store, x(tSpan(1):0.01:tSpan(end))];
    y_store = [y_store, y(tSpan(1):0.01:tSpan(end))];
  
    % re-initialize to the true experimental values
    v0 = V(tSpan(end));
    phi0 = Yaw(tSpan(end));
    x0 = X(tSpan(end));
    y0 = Y(tSpan(end));
end


figure
plot(X_out.data, Y_out.data,'LineWidth',5)
hold on
scatter(x_store, y_store)
legend('Measured','Simulated')
title('Validating kinematic bicycle model')

end