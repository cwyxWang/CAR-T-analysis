function pass = step1_entry(config)
    
    temp = dir(fullfile(config.view_path,'K_*'));
    K_name = {temp.name};
    k_num = size(K_name,2);
    temp = dir(fullfile(config.view_path,'*.dcimg'));
    H_name = {temp.name};
    if size(H_name, 2) == 0
        temp = dir(fullfile(config.view_path,'H_*'));
        H_name = {temp.name};
    end
    h_num = size(H_name,2);
    
    data_num = max(k_num, h_num);
    
    config.timepoint = 0;
    
    for i = 1:data_num
        if i <= h_num
            config.h_path = fullfile(config.view_path, H_name{i});
            fprintf([strrep(config.view_path,'\','/'),'/',H_name{i},'\n']);
        else
            config.h_path = '';
        end
        if i <= k_num
            config.k_path = fullfile(config.view_path, K_name{i});
            fprintf([strrep(config.view_path,'\','/'),'/',K_name{i},'\n']);
        else
            config.k_path = '';
        end
        
        timepoint = step2_thread(config);
        
        config.timepoint = config.timepoint + timepoint;
        if config.actual_pass<=timepoint
            config.actual_pass = 0;
        else
            config.actual_pass = config.pass-timepoint;
        end
    end
    
    pass = config.actual_pass;
end

