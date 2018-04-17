function clist = get_clist(api)
    url =[api, 'local_computer/'];
    %get computer list
    lc_list = webread(url);
    clist.id=[lc_list.objects.id];
    clist.name={lc_list.objects.name};
end