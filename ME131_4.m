% Lab 4
% Front Left encoder works best, then front right, and rears are broken
%% load data from five velocity identification experiments
close all
clear all
clc
exp1 = load('4_2_1650_5.mat', 'sig');
exp2 = load('4_2_1650_4.mat', 'sig');
exp2.sig{16} = exp2.sig{12};
exp3 = load('4_2_1650_8.mat', 'sig');
exp4 = load('4_2_1650_9.mat', 'sig');
exp5 = load('4_2_1650_10.mat', 'sig');
expData = [exp1, exp2, exp3, exp4, exp5];
n=length(expData);

%% wheel radius
radius = 0.05;

%% take the average data from encoders of two front wheels for each experiment
encoderData = {n};
for i = 1:n
    encoderData{i} = (expData(i).sig{1,16}.Data);
end

timeEncoder = {n};
for i = 1:n
    timeEncoder{i} = expData(i).sig{1,16}.Time;
end

%% load the motor pwm data
motorpwmData = {n};
timeMotor = {n};
h = subplot(2,1,1);
hold on
timeSeqMotor = {n};
for i = 1:n
    motorpwmData{i} = expData(i).sig{1,10}.Data - 1500;
    timeMotor{i} = expData(i).sig{1,10}.Time;
    % transfer the unit from day to second
    timeSeqMotor{i} = (timeMotor{i} - timeEncoder{i}(1))*24*3600;
    plot(timeSeqMotor{i}, motorpwmData{i})
    legendInfo{i} = ['test' num2str(i)];
end
legend(legendInfo,'location','best')
xlabel('time [s]')
ylabel('motor pwm')
grid
hold off

%% plot the angular velocity vs. time

subplot(2,1,2)
hold on
angularVelSeq = {n};
timeSeqAngVel = {n};
speedSeq = {n};
for experIndex = 1:n
    timeN = [];
    dataN = [];
    for i = (1+4):4:(length(timeEncoder{experIndex})-4)
        timeN = [timeN, timeEncoder{experIndex}(i)];
        dataN = [dataN, mean(encoderData{experIndex}(i-4:i+4))];
    end
    
    % Estimating the derivative using least-squares polynomial fits
    timeN = (timeN - timeN(1))*24*3600; % unit: second
    dataN = dataN/8*2*pi; % unit: radian
    angularVelocityEst = zeros(1,numel(timeN));
    for i=(1+2):(numel(timeN)-2)
        % Take 5 points, centered at timeN(i), fit with quadratic
        P = polyfit(timeN(i-2:i+2),dataN(i-2:i+2),2);
        % Take derivative of quadratic
        W = polyder(P);
        % Evaluate derivative at timeN(i), to get estimate of angular
        % velocity
        angularVelocityEst(i) = polyval(W, timeN(i));
    end
    angularVelSeq{experIndex} = angularVelocityEst;
    timeSeqAngVel{experIndex} = timeN;
    speedSeq{experIndex} = angularVelSeq{experIndex} * radius;
    plot(timeSeqAngVel{experIndex}, speedSeq{experIndex});
    legendInfo2{experIndex} = ['test' num2str(experIndex)];
end
xlabel('time: [s]')
ylabel('linear velocity: [meter/s]')
grid minor
xlim(get(h, 'XLim'))

%% system ID
vdot = [];
vu = [];
motorpwmData_1 = motorpwmData{1};
motorpwmData_2 = motorpwmData{2};
motorpwmData_3 = motorpwmData{3};
motorpwmData_4 = motorpwmData{4};
motorpwmData_5 = motorpwmData{5};
motorpwmData_1_unique = motorpwmData_1 + linspace(0,1,length(motorpwmData_1))*1E-1;
motorpwmData_2_unique = motorpwmData_2 + linspace(0,1,length(motorpwmData_2))*1E-1;
motorpwmData_3_unique = motorpwmData_3 + linspace(0,1,length(motorpwmData_3))*1E-1;
motorpwmData_4_unique = motorpwmData_4 + linspace(0,1,length(motorpwmData_4))*1E-1;
motorpwmData_5_unique = motorpwmData_5 + linspace(0,1,length(motorpwmData_5))*1E-1;
newmotorpwmData = {diag(motorpwmData_1_unique), diag(motorpwmData_2_unique), diag(motorpwmData_3_unique), diag(motorpwmData_4_unique), diag(motorpwmData_5_unique)};

timeSeqMotor_1 = timeSeqMotor{1};
timeSeqMotor_2 = timeSeqMotor{2};
timeSeqMotor_3 = timeSeqMotor{3};
timeSeqMotor_4 = timeSeqMotor{4};
timeSeqMotor_5 = timeSeqMotor{5};
timeSeqMotor_1_unique = timeSeqMotor_1 + linspace(0,1,length(timeSeqMotor_1))*1E-1;
timeSeqMotor_2_unique = timeSeqMotor_2 + linspace(0,2,length(timeSeqMotor_2))*1E-1;
timeSeqMotor_3_unique = timeSeqMotor_3 + linspace(0,3,length(timeSeqMotor_3))*1E-1;
timeSeqMotor_4_unique = timeSeqMotor_4 + linspace(0,4,length(timeSeqMotor_4))*1E-1;
timeSeqMotor_5_unique = timeSeqMotor_5 + linspace(0,5,length(timeSeqMotor_5))*1E-1;
newtimeSeqMotor = {diag(timeSeqMotor_1_unique), diag(timeSeqMotor_2_unique), diag(timeSeqMotor_3_unique), diag(timeSeqMotor_4_unique), diag(timeSeqMotor_5_unique)};
for i = 1:n
    newmotorpwmData{i} = interp1(newtimeSeqMotor{i}, newmotorpwmData{i}, timeSeqAngVel{i},...
    'linear','extrap');
    startingIndex{i} = find(speedSeq{i}>0,1);
    vdot = [vdot; ...
    (diff(speedSeq{i}(startingIndex{i}:end))./diff(timeSeqAngVel{i}(startingIndex{i}:end)))'];
    vu = [vu; speedSeq{i}(startingIndex{i}:end-1)', newmotorpwmData{i}(startingIndex{i}:end-1)'];
end
% Do least square to get the transfer function:V(s)/U(s) = b/(s ? a).
ab = vu\vdot;
a = ab(1);
b = ab(2);
model = tf([b], [1 -a])

%% plot the simulation of the identified model
u =@(t) interp1(timeSeqAngVel{1}, newmotorpwmData{1}, t);
[tODE, vODE] = ode45(@(t,v) dyn(t, v, u, ab), ...
[timeSeqAngVel{1}(startingIndex{1}),timeSeqAngVel{1}(end)], 0);
plot(tODE, vODE, 'linewidth',2)
legendInfo2{n+1} = ['Fitted system'];
legend(legendInfo2,'location', 'best')
hold off

%%-------------------------------------------------------------------------
%% load data from five velocity identification experiments
close all
clear all
clc
exp1 = load('braking_1440_1.mat', 'sig');
% exp1.sig{16} = exp1.sig{12};
exp2 = load('braking_1440_2.mat', 'sig');
% exp2.sig{16} = exp2.sig{12};
exp3 = load('braking_1440_3.mat', 'sig');
% exp3.sig{16} = exp3.sig{12};
exp4 = load('braking_1440_4.mat', 'sig');
exp5 = load('brake_1440_5.mat', 'sig');
% exp5.sig{16} = exp5.sig{12};
expData = [exp1, exp2, exp3, exp4, exp5];
n=length(expData);

%% wheel radius
radius = 0.05;

%% take the average data from encoders of two front wheels for each experiment
encoderData = {n};
for i = 1:n
    encoderData{i} = (expData(i).sig{1,16}.Data);
end

timeEncoder = {n};
for i = 1:n
    timeEncoder{i} = expData(i).sig{1,12}.Time;
end

%% load the motor pwm data
motorpwmData = {n};
timeMotor = {n};
h = subplot(2,1,1);
hold on
timeSeqMotor = {n};
for i = 1:n
    motorpwmData{i} = expData(i).sig{1,10}.Data - 1500;
    timeMotor{i} = expData(i).sig{1,10}.Time;
    % transfer the unit from day to second
    timeSeqMotor{i} = (timeMotor{i} - timeEncoder{i}(1))*24*3600;
    plot(timeSeqMotor{i}, motorpwmData{i})
    legendInfo{i} = ['test' num2str(i)];
end
legend(legendInfo,'location','best')
xlabel('time [s]')
ylabel('motor pwm')
grid
hold off

%% plot the angular velocity vs. time

subplot(2,1,2)
hold on
angularVelSeq = {n};
timeSeqAngVel = {n};
speedSeq = {n};
for experIndex = 1:n
    timeN = [];
    dataN = [];
    for i = (1+4):4:(length(timeEncoder{experIndex})-4)
        timeN = [timeN, timeEncoder{experIndex}(i)];
        dataN = [dataN, mean(encoderData{experIndex}(i-4:i+4))];
    end
    
    % Estimating the derivative using least-squares polynomial fits
    timeN = (timeN - timeN(1))*24*3600; % unit: second
    dataN = dataN/8*2*pi; % unit: radian
    angularVelocityEst = zeros(1,numel(timeN));
    for i=(1+2):(numel(timeN)-2)
        % Take 5 points, centered at timeN(i), fit with quadratic
        P = polyfit(timeN(i-2:i+2),dataN(i-2:i+2),2);
        % Take derivative of quadratic
        W = polyder(P);
        % Evaluate derivative at timeN(i), to get estimate of angular
        % velocity
        angularVelocityEst(i) = polyval(W, timeN(i));
    end
    angularVelSeq{experIndex} = angularVelocityEst;
    timeSeqAngVel{experIndex} = timeN;
    speedSeq{experIndex} = angularVelSeq{experIndex} * radius;
    plot(timeSeqAngVel{experIndex}, speedSeq{experIndex});
    legendInfo2{experIndex} = ['test' num2str(experIndex)];
end
xlabel('time: [s]')
ylabel('linear velocity: [meter/s]')
grid minor
xlim(get(h, 'XLim'))

%% system ID
vdot = [];
vu = [];
motorpwmData_1 = motorpwmData{1};
motorpwmData_2 = motorpwmData{2};
motorpwmData_3 = motorpwmData{3};
motorpwmData_4 = motorpwmData{4};
motorpwmData_5 = motorpwmData{5};
motorpwmData_1_unique = motorpwmData_1 + linspace(0,1,length(motorpwmData_1))*1E-1;
motorpwmData_2_unique = motorpwmData_2 + linspace(0,1,length(motorpwmData_2))*1E-1;
motorpwmData_3_unique = motorpwmData_3 + linspace(0,1,length(motorpwmData_3))*1E-1;
motorpwmData_4_unique = motorpwmData_4 + linspace(0,1,length(motorpwmData_4))*1E-1;
motorpwmData_5_unique = motorpwmData_5 + linspace(0,1,length(motorpwmData_5))*1E-1;
newmotorpwmData = {diag(motorpwmData_1_unique), diag(motorpwmData_2_unique), diag(motorpwmData_3_unique), diag(motorpwmData_4_unique), diag(motorpwmData_5_unique)};

timeSeqMotor_1 = timeSeqMotor{1};
timeSeqMotor_2 = timeSeqMotor{2};
timeSeqMotor_3 = timeSeqMotor{3};
timeSeqMotor_4 = timeSeqMotor{4};
timeSeqMotor_5 = timeSeqMotor{5};
timeSeqMotor_1_unique = timeSeqMotor_1 + linspace(0,1,length(timeSeqMotor_1))*1E-1;
timeSeqMotor_2_unique = timeSeqMotor_2 + linspace(0,2,length(timeSeqMotor_2))*1E-1;
timeSeqMotor_3_unique = timeSeqMotor_3 + linspace(0,3,length(timeSeqMotor_3))*1E-1;
timeSeqMotor_4_unique = timeSeqMotor_4 + linspace(0,4,length(timeSeqMotor_4))*1E-1;
timeSeqMotor_5_unique = timeSeqMotor_5 + linspace(0,5,length(timeSeqMotor_5))*1E-1;
newtimeSeqMotor = {diag(timeSeqMotor_1_unique), diag(timeSeqMotor_2_unique), diag(timeSeqMotor_3_unique), diag(timeSeqMotor_4_unique), diag(timeSeqMotor_5_unique)};
for i = 1:n
    newmotorpwmData{i} = interp1(newtimeSeqMotor{i}, newmotorpwmData{i}, timeSeqAngVel{i}, 'linear','extrap');
    brakingIndex{i} = find(newmotorpwmData{i}<0,1);
    vdot = [vdot; ...
    (diff(speedSeq{i}(brakingIndex{i}:end))./diff(timeSeqAngVel{i}(brakingIndex{i}:end)))'];
    vu = [vu; speedSeq{i}(brakingIndex{i}:end-1)', newmotorpwmData{i}(brakingIndex{i}:end-1)'];
end
% Do least square to get the transfer function:V(s)/U(s) = b/(s ? a).
ab = vu\vdot;
a = ab(1);
b = ab(2);
model = tf([b], [1 -a])

%% plot the simulation of the identified model
u =@(t) interp1(timeSeqAngVel{1}, newmotorpwmData{1}, t);
[tODE, vODE] = ode45(@(t,v) dyn(t, v, u, ab), ...
[timeSeqAngVel{1}(brakingIndex{1}),timeSeqAngVel{1}(end)], speedSeq{1}(brakingIndex{1}));
plot(tODE, vODE, 'linewidth',2)
legendInfo2{n+1} = ['Fitted system'];
legend(legendInfo2,'location', 'best')
hold off

%%-------------------------------------------------------------------------
%% load data from steering mapping experiments
clear all
clc
dinfo = dir('steering/*.mat'); % just get all steering data files
for K = 1 : length(dinfo)
    thisfile = dinfo(K).name
    expData{K}=load(strcat('steering/',thisfile),'sig');
end
clear dinfo
n=K;
clear K

%% parameters
radius = 0.05;
wheelbase = 0.3;

%% load time data for five experiments
% take the average data from encoders of two front wheels for each
% experiment.
encoderData = {n};
timeEncoder = {n};
% Change this for our own encoder data slots
for i = 1:n
    encoderData{i} = (expData{i}.sig{1,12}.Data +...
        expData{i}.sig{1,13}.Data)/2;
    timeEncoder{i} = expData{i}.sig{1,12}.Time;
end

%% load yaw for five experiments
timeYaw = {n};
for i = 1:n
    timeYaw{i} = expData{i}.sig{1,3}.Time;
end

% take the average data from encoders of two front wheels for each
% experiment.
yawData = {n};
for i = 1:n
    yawData{i} = unwrap(expData{i}.sig{1,3}.Data);
end

%% load the servo pwm data
servopwmData = {n};
timeServo = {n};
h = subplot(3,1,1);
hold on
timeSeqServo = {n};
for i = 1:n
    servopwmData{i} = expData{i}.sig{1,11}.Data;
    timeServo{i} = expData{i}.sig{1,11}.Time;
    % transfer the unit from day to second
    timeSeqServo{i} = (timeServo{i} - timeEncoder{i}(1))*24*3600;
    plot(timeSeqServo{i}, servopwmData{i})
    five_index = find(timeSeqServo{i}>2.5,1);
    servo_input{i} = servopwmData{i}(five_index);
end

xlabel('time [s]')
ylabel('servo pwm')
xlim([0 6])
grid
hold off

%% plot the angular velocity vs. time
subplot(3,1,2)
hold on
angularVelSeq = {n};
timeSeqAngVel = {n};
speedSeq = {n};
speed = {n};
for experIndex = 1:n
    timeN = [];
    dataN = [];
    for i = (1+4):4:(length(timeEncoder{experIndex})-4)
        timeN = [timeN, timeEncoder{experIndex}(i)];
        dataN = [dataN, mean(encoderData{experIndex}(i-4:i+4))];
    end
    % Estimating the derivative using least-squares polynomial fits
    timeN = (timeN - timeN(1))*24*3600; % unit: second
    dataN = dataN/8*2*pi; % unit: radian
    angularVelocityEst = zeros(1,numel(timeN));
    for i=(1+2):(numel(timeN)-2)
        % Take 5 points, centered at timeN(i), fit with quadratic
        P = polyfit(timeN(i-2:i+2), dataN(i-2:i+2),2);
        % Take derivative of quadratic
        W = polyder(P);
        % Evaluate derivative at timeN(i), to get estimate of angular
        % velocity
        angularVelocityEst(i) = polyval(W,timeN(i));
    end
    angularVelSeq{experIndex} = angularVelocityEst;
    timeSeqAngVel{experIndex} = timeN;
    speedSeq{experIndex} = angularVelSeq{experIndex}*radius;
    
    five_index = find(timeSeqAngVel{experIndex}>2.5,1);
    speed{experIndex} = speedSeq{experIndex}(five_index);
    plot(timeSeqAngVel{experIndex}, speedSeq{experIndex});
end
xlabel('time: [s]')
ylabel('linear velocity: [meter/s]')
grid minor

%% Get yaw rate data
subplot(3,1,3)
hold on
yawSeq = {n}
yawRate = {n};
timeSeqYawRate = {n};
for experIndex = 1:n
    timeN = [];
    dataN = [];
    for i = (1+4):4:(length(timeYaw{experIndex})-4)
        timeN = [timeN, timeYaw{experIndex}(i)];
        dataN = [dataN, mean(yawData{experIndex}(i-4:i+4))];
    end
    % Estimating the derivative using least-squares poly fits
    timeN = (timeN - timeN(1))*24*3600; % unit: second
    timeSeqYawRate{experIndex} = timeN;
    % yawRate as slope at 5 seconds
    yawSeq{experIndex}=dataN;
    five_indx = find(timeSeqYawRate{experIndex}>2.5,1);
    yawRate{experIndex} = (yawSeq{experIndex}(five_index)-yawSeq{experIndex}(five_index-1))/...
        (timeSeqYawRate{experIndex}(five_index)-timeSeqYawRate{experIndex}(five_index-1));
    plot(timeSeqYawRate{experIndex},yawSeq{experIndex});
end
xlabel=('time: [s]')
ylabel('linear velocity: [meter/s]')
grid minor
xlim(get(h, 'Xlim'))

%% Convert velocity and yaw information for each sample
% steering angle = atan(L*yawrate/vel)
steering_angle={n};
figure
hold on
for experIndex=1:n
    steering_angle{experIndex} = wheelbase * yawRate{experIndex} / speed{experIndex};
    if servo_input{experIndex}>1000
        if experIndex~=13
            scatter(servo_input{experIndex}, steering_angle{experIndex}, 'kx')
        end
    end
end
xlabel('servo input PWM')
ylabel('steering angle')
title('servo pwm -> steering angle mapping')

%% plot the line of best fit
x1=[];
y1=[];
for i=1:n
    x1 = [x1 servo_input{i}];
    y1 = [y1 steering_angle{i}];
end
P = polyfit(x1, y1, 1);
yfit1 = P(1)*x1+P(2);
plot(x1, yfit1, 'b-');


%% 4.1
% Run thecheckMovementDemoOutput.mscript to validate your recorded data
% from each experiment.  Submit a screenshot of the command window with the console output.

%% 4.2 - Velocity
% Find the transfer function in the frequency domain,Gpwm?v(s).
% A plot comparing your first-order model and your five recorded velocity vs.  time datasets.

% car code: 1120632. experiment 32, 33, 34, and 35 [4.2.1, 4.2.2, 4.2.4,
% 4.2.5]. You must manually select these experiments in order as you
% interact with the download_from_cloud script. put 0 [no] and 0 [no] for
% plotting and viewing video for each experiment. There are only four
% trials b/c 4.2.3 didn't record properly, we may want to redo it.
C = {};
for i = 49:54
    C{end+1} = exp{i};
end

req_data = {};
counter = 0;
for j = 1:5
    vel_x = C{j}{12};
    acc_long = C{j}{7};
    servo_pwm = C{j}{11};
    motor_pwm = C{j}{10};
    
    acc_long_rs = resample(acc_long, vel_x.time, 'linear');
    servo_pwm_rs = resample(servo_pwm, vel_x.time, 'linear');
    motor_pwm_rs = resample(motor_pwm, vel_x.time, 'linear');
    
    motor_pwm_norm = motor_pwm_rs.data - 1500;
    
    v_matrix = [vel_x.data, motor_pwm_norm];
    t = (C{j}{12}.Time - C{j}{12}.Time(1))*86400;
    
    for k = 1:length(v_matrix)
        if isnan(v_matrix(k,2))
            v_matrix = v_matrix(1:k-1,:);
            acc_long_rs_data = acc_long_rs.data(1:k-1,:);
            servo_pwm_rs_data = servo_pwm_rs.data(1:k-1,:);
            motor_pwm_norm = motor_pwm_norm(1:k-1,:);
            t = t(1:k-1,:);
            break;
        end
    end
%     whos v_matrix
%     whos acc_long_rs_data
    
    if length(v_matrix) ~= length(acc_long_rs_data) & length(v_matrix) > length(acc_long_rs_data)
        v_matrix = v_matrix(1:length(acc_long_rs_data),:);
        servo_pwm_rs_data = servo_pwm_rs.data(1:length(acc_long_rs_data),:);
        motor_pwm_norm = motor_pwm_norm(1:length(acc_long_rs_data),:);
        t = t(1:length(acc_long_rs_data),:);
    end
    if length(v_matrix) ~= length(acc_long_rs_data) & length(v_matrix) < length(acc_long_rs_data)
        acc_long_rs_data = acc_long_rs_data(1:length(v_matrix),:);
        servo_pwm_rs_data = servo_pwm_rs.data(1:length(v_matrix),:);
        motor_pwm_norm = motor_pwm_norm(1:length(v_matrix),:);
        t = t(1:length(v_matrix),:);
    end
        
    x = v_matrix\acc_long_rs_data;
    
    a = x(1);
    b = x(2);
    
    model = tf([b], [1 -a]);
%     t = (C{j}{12}.Time - C{j}{12}.Time(1))*86400;
    
    [y,t,x] = lsim(model,motor_pwm_norm,linspace(t(1),t(end),length(t)));
    
    indiv_arr = {};
    indiv_arr{end+1} = vel_x;
    indiv_arr{end+1} = acc_long_rs_data;
    indiv_arr{end+1} = servo_pwm_rs_data;
    indiv_arr{end+1} = motor_pwm_norm;
    indiv_arr{end+1} = t;
    indiv_arr{end+1} = y;
    
    req_data{end+1} = indiv_arr;
    counter = counter + 1;
end
%%
% Some stuff that was being tested
    vel_x_32 = exp{32}{12};
    vel_x_33 = exp{33}{12};
    vel_x_34 = exp{34}{12};
    vel_x_35 = exp{35}{12};

    % Exp 34 yields best fit line
    exp_num=34;
    
    vel_x = exp{exp_num}{12};
    acc_long = exp{exp_num}{7};
    servo_pwm = exp{exp_num}{11};
    motor_pwm = exp{exp_num}{10};
    
    acc_long_rsmpld = resample(acc_long, vel_x.time, 'linear');
    servo_pwm_rsmpld = resample(servo_pwm, vel_x.time, 'linear');
    motor_pwm_rsmpld = resample(motor_pwm, vel_x.time, 'linear');
    
    v_matrix = [vel_x.data, motor_pwm_rsmpld.data-1500];
    
    for z = 1:length(v_matrix)
        if isnan(v_matrix(z,2))
            v_matrix = v_matrix(1:z-1,:);
            acc_long_rsmpld = acc_long_rsmpld.data(1:z-1,:);
            motor_pwm_rsmpld = motor_pwm_rsmpld.data(1:z-1,:);
            break;
        end
    end
    
    x = v_matrix\acc_long_rsmpld;
    
    a = x(1);
    b = x(2);
    
    model = tf([b], [1, -a]);

    t_main = (exp{exp_num}{12}.Time - exp{exp_num}{12}.Time(1))*86400;
    t = t_main(1:length(motor_pwm_rsmpld));
    
    % Deliverable 2
    figure(1)
    plot(t(1:length(motor_pwm_rsmpld)),motor_pwm_rsmpld)
    
    % Determine which is the best fit line to use
%     figure(2)
    [y,t,x] = lsim(model,motor_pwm_rsmpld-1500,linspace(t(1),t(end),length(t)));
%     plot(t,y)
    best_fit_line = y;
%     plot(t,best_fit_line)
    % Use exp_num = 34
    
%     figure(3)
%     plot(movmean(vel_x.data,10))

figure(2)
hold on
plot(best_fit_line)
plot(movmean(vel_x_32.data,30))
plot(movmean(vel_x_33.data,30))
plot(movmean(vel_x_34.data,30))
plot(movmean(vel_x_35.data,30))


% for exp_num = 50:54
%     vel_x = exp{exp_num}{12};
%     acc_long = exp{exp_num}{7};
%     servo_pwm = exp{exp_num}{11};
%     motor_pwm = exp{exp_num}{10};
%     
%     acc_long_rsmpld = resample(acc_long, vel_x.time, 'linear');
%     servo_pwm_rsmpld = resample(servo_pwm, vel_x.time, 'linear');
%     motor_pwm_rsmpld = resample(motor_pwm, vel_x.time, 'linear');
%     
%     v_matrix = [vel_x.data, motor_pwm_rsmpld.data-1500];
%     
%     for z = 1:length(v_matrix)
%         if isnan(v_matrix(z,2))
%             acc_long_rsmpld = acc_long_rsmpld.data(1:z-1,:);
%             motor_pwm_rsmpld = motor_pwm_rsmpld.data(1:z-1,:);
%             break;
%         end
%     end
%     
%     x = v_matrix\acc_long_rsmpld;
%     
%     a = x(1);
%     b = x(2);
%     
%     model = tf([b], [1, -a]);
% 
%     t = (exp{exp_num}{12}.Time - exp{exp_num}{12}.Time(1))*86400;
%     t = t(1:length(motor_pwm_rsmpld));
%     
%     % Deliverable 2
%     figure(1)
%     plot(t(1:length(motor_pwm_rsmpld)),motor_pwm_rsmpld)
%     
%     % Deliverable 3
%     figure(2)
%     [y,t,x] = lsim(model,motor_pwm_rsmpld-1500,linspace(t(1),t(end),length(t)));
%     plot(t,y)
%     
%      
% vel_x = exp{50}{12};
% acc_long = exp{50}{7};
% servo_pwm = exp{50}{11};
% motor_pwm = exp{50}{10};
% 
% acc_long_rsmpld = resample(acc_long, vel_x.time, 'linear');
% servo_pwm_rsmpld = resample(servo_pwm, vel_x.time, 'linear');
% motor_pwm_rsmpld = resample(motor_pwm, vel_x.time, 'linear');
% 
% v_matrix = [vel_x.data, motor_pwm_rsmpld.data-1500];
% length(v_matrix)
% 
% for z=1:191
%     if isnan(v_matrix(z,2))
%         v_matrix = v_matrix(1:z-1,:);
%         acc_long_rsmpld = acc_long_rsmpld.data(1:z-1,:);
%         motor_pwm_rsmpld = motor_pwm_rsmpld.data(1:z-1,:);
%         break;
%     end
% end
% 
% 
% 
% x = v_matrix\acc_long_rsmpld;
% 
% a = x(1)
% b = x(2)
% 
% model = tf([b], [1, -a]);
% 
% t = (exp{50}{12}.Time - exp{50}{12}.Time(1))*86400;
% t = t(1:length(motor_pwm_rsmpld));
% 
% figure(1)
% plot(t(1:length(motor_pwm_rsmpld)),motor_pwm_rsmpld)
% 
% figure(2)
% [y,t,x] = lsim(model,motor_pwm_rsmpld-1500,linspace(t(1),t(end),length(t)));
% plot(t,y)

% clear all; close all;
% 
% 
% 
% clear all; close all;
% a = []; v = [];
% for idx = 1:1% 4
%     download_from_cloud 
%     a_raw = sig{7};
%     a_smooth = timeseries(movmean(a_raw.Data,10), a_raw.time); % smooth accelerometer data
%     u_raw = sig{10};
%     %v_raw = sig{12};
%     %v_smooth = timeseries(movmean(v_raw.Data,10), v_raw.time); % smooth velocity
%     
%     enc_resamp = resample(sig{16}.data, time, 'linear');
%     ang_vel = [];
%     for idx = 1:numel(time) - 1
%         enc_diff = enc_resamp(idx+1) - enc_resamp(idx);
%         time_diff = time(idx+1) - time(idx);
%         ang_vel = [ang_vel enc_diff/time_diff];
%     end
%     r = 5; % 5cm radius
%     vel_raw = r .* ang_vel / 8;
%     vel_raw(numel(time)) = 0;
%     %plot(time, movmean(vel, 10))
%     vel = movmean(vel_raw, 10);
%     vel_ts = timeseries(vel', sig{12}.time); %sig{12}.time or time?
% 
%     u_resample = resample(u_raw, sig{12}.time, 'linear');
%     a_resample = resample(a_smooth, sig{12}.time, 'linear');
%     %v_resample = v_raw;
%     %find the first NaN in each column
%     u_NaN = find(isnan(u_resample.data), 1); 
%     a_NaN = find(isnan(a_resample.data), 1);
%     %v_NaN = find(isnan(v_resample.data), 1);
%     Trun = min([u_NaN a_NaN])
%     % truncate all columns at the first NaN
%     u_seg = timeseries(u_resample.data(1:Trun-1), u_resample.time(1:Trun-1))
%     a_seg = timeseries(a_resample.data(1:Trun-1), a_resample.time(1:Trun-1));
%     v_seg = timeseries(vel_ts.data(1:Trun-1), vel_ts.time(1:Trun-1))
%     u_seg.data = u_seg.data - 1500; % 1500 = 0 motor torque
%     for kdx = 1:numel(u_seg.data)
%         if u_seg.data(kdx) < 0
%             u_seg.data(kdx) = 0
%         end
%     end
%     a = [a; a_seg.data];
%     %find(isnan(a), 1); 
%     v = [v; v_seg.data u_seg.data];
%     %find(isnan(v), 1); 
% 
% end

%[a;b] = v\a

%a =  -0.007941771812576, b = 0.000017380723421
% a = 0.013452897943938, b = 0.000000142533124  with self-derived velocity
% a = -0.0342 b = 0.0078
% a = 0.078, b = 0.342 guessed, gives better velocity profile
% A plot of motor actuation PWM vs.  time data
%plot(u_raw)


%% 4.3 - Braking



