% Functions and script for processing data from Dator.
% Set api url
api = 'http://dator.forge9.com/api/v1/';


% loc_comp_id = input('Type your 7 digit computer ID (read it from cloud.cfg file on your BARC)?');
loc_comp_id = 1120632;

% get experiment list
elist = get_explist(api,loc_comp_id);
disp('----EXPERIMENT LIST----')
for k=1:length(elist.name)
    disp([num2str(k),' - ',elist.name{k}])
end
% select experiment here 
exp_idx = input('Which experiment you want to download?');
disp('...downloading..be patient...')
exp_name=elist.name{exp_idx};
exp_id=elist.id(exp_idx);


% get experiment signals
sig = get_signal(loc_comp_id,exp_name);

%% remove multiple experiments from signal
newsig = sig;
for i = 1:length(sig)
    time = sig{1, i}.Time;
    data = sig{1, i}.Data;
    time = time*3600*24;
    time= time-time(1);
    TS = time(2:end)-time(1:end-1);
    bigdelay = find(TS>1); %sample rate should not be greater than 1Hz
    if ~isempty(bigdelay)
        lastdelay = bigdelay(end)+1;
        ts= timeseries(sig{1, i}.Data(lastdelay:end),sig{1, i}.Time(lastdelay:end),'Name',sig{1, i}.Name);
        newsig{1,i} = ts;
    end
end
sig = newsig;
%%
% get experiment media
exp_media = get_exp_media(exp_id);
savevideo=1;
if savevideo
    if ~isempty(exp_media)
        outfilename = websave('local_avi.avi',exp_media);
        vidObj=VideoReader('local_avi.avi');
    else
        vidObj=[];
    end
end

% plot experiment signals
close all
yanswer = input('Want to plot all signals? (1-Yes, 0-No)');
if yanswer
    for j=1:length(sig)
        s=sig{j};
        figure
        plot(s)
    end
end


%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%play video
yanswer = input('Want to play video ? (1-Yes, 0-No)');
if yanswer
    if ~isempty(exp_media)
        vidObj=VideoReader('local_avi.avi');
        figure(1)
        ax=subplot(2,1,1);
        currAxes = ax;
        bx=subplot(2,1,2);
        hold on
        s=sig{1};
        plot(s)
        tv=s.Time(1);
        while hasFrame(vidObj)
            vidFrame = readFrame(vidObj);
            image(vidFrame, 'Parent', currAxes);
            %currAxes.Visible = 'off';
            tv=datenum(tv+seconds(1/vidObj.FrameRate));
            ts1 = resample(s,tv);
            axes(bx)
            plot(ts1,'ro')
            ylabel(bx,s.Name)
            pause(1/vidObj.FrameRate);
        end
    end
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
%% function list

function clist = get_clist(api)
    url =[api, 'local_computer/'];
    %get computer list
    lc_list = webread(url);
    clist.id=[lc_list.objects.id];
    clist.name={lc_list.objects.name};
end

function elist = get_explist(api,lc_id)
    % get experiment list for a computer
    url = [api,'experiment/?local_computer_id=',num2str(lc_id)];
    exp_data = webread(url);
    elist.id=[exp_data.objects.id];
    elist.name={exp_data.objects.name};
end

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

function exp_media = get_exp_media(exp_id)
    url = ['http://dator.forge9.com/api/v1/setting/?experiment_id=',num2str(exp_id)];
    out = webread(url);
    try
        exp_media=out.objects.value;
    catch
        warning('No video associted with this experiment')
        exp_media=[];
    end
end