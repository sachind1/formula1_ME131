function elist = get_explist(api,lc_id)
    % get experiment list for a computer
    url = [api,'experiment/?local_computer_id=',num2str(lc_id)];
    exp_data = webread(url);
    elist.id=[exp_data.objects.id];
    elist.name={exp_data.objects.name};
end
