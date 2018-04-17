clc
clear
close all
download_from_cloud

% Tony's BARC ID: 1120625
% Check experiment 12 (Movement_Demo_01) or 13(Movement_Demo_02) for the
% an example of what the test results should look like

correct_sensors = {'roll';
                    'pitch';
                    'yaw';
                    'angular_velocity_x';
                    'angular_velocity_y';
                    'angular_velocity_z';
                    'linear_acceleration_x';
                    'linear_acceleration_y';
                    'linear_acceleration_z';
                    'motor_pwm';
                    'servo_pwm';
                    'encoder_FL';
                    'encoder_FR';
                    'encoder_BL';
                    'encoder_BR';
                    'velocity_FL';
                    'velocity_FR';
                    'velocity_BL';
                    'velocity_BR';
};
%% Creates a cell array with the names in the correct order.
% If the experiment name is used for multiple runs, data is appended to
% that sig cell array. This will remove everything before the last trial
sig = fix_Multiple_Experiment_Sig(sig);

user_sensors = cell(1,length(sig));
for i = 1:length(sig)
    user_sensors{i} = sig{i}.Name;
end
correct_sig = cell(1,length(correct_sensors));
for i = 1:length(correct_sensors)
    if ~isempty(find(contains(user_sensors,correct_sensors{i})))
        correct_sig{i} = [sig{find(contains(user_sensors,correct_sensors{i}))}];
    else
        correct_sig{i} = timeseries(false,1,'Name',['Error: MISSING ' correct_sensors{i}]);
        
    end
end         
sig = correct_sig; 
%%
for j=1:length(sig)
    try
        if ~isequal(class(sig{j}.Data),'logical')
            if j<4
                sig{j}.Data = unwrap(sig{j}.Data);
            end

            if isequal(sig{j}.Name,'encoder_FL')
                encoderStart = j;
            end
            correct_sensors(find(contains(correct_sensors,sig{j}.Name))) = [];
        end
    catch
        disp('Signal cell array is not ok');
    end
end

if isempty(correct_sensors)
    disp('All sensors are enabled.')
    correct = true;
else
    disp('Missing the following sensor outputs')
    for i=1:length(correct_sensors)
        disp(correct_sensors{i});
    end
    disp('--------------------------------')
    correct = false;
end

if correct == true
    for i = 1:length(sig)
        % Check if signal is empty
        if isempty(sig{i}.Data)
            fprintf([sig{i}.Name,' is empty.\n'])
            correct = false;
        end

        % Check if signal is the correct length
        if length(sig{i}.Data)<300
            fprintf([sig{i}.Name,' should have more data points.\n'])
            correct = false;
        end
    end

    % Check that encoders are working
    if i>=encoderStart && i<=encoderStart+3 && (mean(sig{i}.Data)< 20)
        fprintf([sig{i}.Name,' not working properly.\n'])
        correct = false;
    end

    if ~isempty(vidObj)
        try
            vidFrame = readFrame(vidObj);
        end
        if mean(mean(vidFrame(:,:,1)))<20 &&mean(mean(vidFrame(:,:,2)))<20 &&mean(mean(vidFrame(:,:,3)))<20
            fprintf('Video black\n')
            correct = false;
        end
    else
        fprintf('Video is missing\n')
        correct = false;
    end

    if max(sig{3}.Data)-min(sig{3}.Data)<1
        yawrange = max(sig{3}.Data)-min(sig{3}.Data);
        fprintf('Not enough yaw from turning. Your yaw range is %0.4f \n',yawrange)
        correct = false;
    end


    if (mean(sig{9}.Data)< 0) && isequal(sig{9}.Name,'linear_acceleration_z')
        % Check that imu is mounted in the correct orientation
        fprintf('IMU mounted upside-down. Value should be approximately +9.81\n')
        correct = false;
    end
    if max(sig{11}.Data)~=1800 || min(sig{11}.Data)~=1200
        fprintf('Servo PWM not correct.Min = 1200 & Max = 1800 but you got Min: %d & Max: %d\n',min(sig{11}.Data), max(sig{11}.Data))
        correct = false;
    end
    
    last_pct = 10;
    lenAccel = length(sig{9}.Data);
    accelFinalRange = max(sig{9}.Data(lenAccel-round(lenAccel/last_pct):end))-min(sig{9}.Data(lenAccel-round(lenAccel/last_pct):end));
    if accelFinalRange>0.4
        fprintf('Final z-acceleration not ok! Wait for process to finish cleanly before moving BARC when recording experiment\n',accelFinalRange)
        correct = false;
    end
end


if correct == true
    disp('Data is correct.')
end