function sig = get_signal(lc_id,exp_name)
    % get signal from an experiment
    %http://dator.forge9.com/data_api/v1/local_computer/1120540/find_signals/?experiment=test_001_2016.04.15_14.13_2016.04.15_14.13&include_data
    url = ['http://dator.forge9.com/data_api/v1/local_computer/',num2str(lc_id),'/find_signals/?experiment=',exp_name,'&include_data'];
    exp_data = webread(url);
    signal_car=[];
    close all
%   basetime=datenum8601(exp_data(j).fields.created_at);
    basetime=datenum(1969,1,1);
    for j=1:length(exp_data)
        % make timevec
        if ~isempty(exp_data(j).data)
            timevec=datenum(basetime+seconds(exp_data(j).data(:,1)));
            s=timeseries(exp_data(j).data(:,2),timevec,'Name',exp_data(j).fields.name);
            sig{j}=s;
        else
            sig{j}=timeseries([],[],'name',exp_data(j).fields.name);
            disp(['warning!! empty signal: ',exp_data(j).fields.name])
        end
    end
end
