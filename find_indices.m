for i = 1:length(expData)
    for j = 1:length(expData{i}.sig)
        if expData{i}.sig{j}.name == 'encoder_FL'
            fprintf('%d %d encoder_FL', i, j)
        elseif expData{i}.sig{j}.name == 'encoder_FR'
            fprintf('%d %d encoder_FR', i, j)
        elseif expData{i}.sig{j}.name == 'yaw'
            fprintf('%d %d yaw', i, j)
        elseif expData{i}.sig{j}.name == 'servo_pwm'
            fprintf('%d %d servo_pwm', i, j)
        else
            fprintf('%d has none of these', i)
        end
    end
end