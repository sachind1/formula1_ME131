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