function configure(config)
% v_12_30_2-6_ds：12个细胞（n_split）；每个体积拍30张（即一次volume_num过程中拍30张，n_slice_per_stack），2-6（2.6，光片扫描步长）
% preview：只跑MIP，只跑投影   mip_num：校准体积数   mean：噪声   std：标准差   flipz：true or false
% wview：滨松拍出来的是否自动转上下画幅   volume_num：（4，一个细胞拍4个来回），gpu_index：选择用哪个GPU（该参数可删除）
% rotate：转平，布尔值，否的话就是要中间过程未转平的结果
    t1 = clock;
    % 读根目录下的所有v_开头的文件名
    temp = dir(fullfile(config.root,'v_*'));
    View_name = {temp.name};
    view_num = size(View_name, 2); 
    pass = config.pass;
    config.roi_num = 1;
    config.pixelsize = config.pixelsize*config.binning;

    for i = 1:view_num
        
        config.pass = pass/config.roi_num;
        config_str = split(View_name{i},'_');
        config.volume_num = str2double(config_str(2));
        config.roi_num = str2double(config_str(3));
        config.actual_test_num = config.test_num*config.roi_num*config.volume_num;
        config.actual_pass = config.pass*config.roi_num;
        config.slice_per_stack = str2double(config_str(4));
        config.stepsize = str2double(strrep(config_str(5),'-','.'));
        config.factor = config.stepsize / config.pixelsize;
        config.view_path = fullfile(config.root,View_name{i});
        
        pass = step1_entry(config);
        make_montage(config.view_path);
    end
    t2 = clock;
    dt1 = etime(t2,t1);
    
    if config.segment
        data_path = config.root;
        ckpt = config.ckpt;
        thread_num = config.gpu_num;
        if thread_num > 1
            parfor i = 0:thread_num-1
                run_3dunet(data_path, ckpt, i, thread_num);
            end
        else
            run_3dunet(data_path, ckpt, 0, thread_num);
        end
        
        t3 = clock;
        dt2 = etime(t3,t2);
    end
    
    if config.analysis
        step1_setup_analysis(config.root, config.volume_num, 0, false);
        t4 = clock;
        dt3 = etime(t4,t3);
    end
    
    m = floor(dt1/60);
    s = dt1-m*60;
    h = floor(m/60);
    m = m-h*60;
    fprintf(['重建总耗时',num2str(h),'h',num2str(m),'m',num2str(s),'s\n'])
    if config.segment
        m = floor(dt2/60);
        s = dt2-m*60;
        h = floor(m/60);
        m = m-h*60;
        fprintf(['分割总耗时',num2str(h),'h',num2str(m),'m',num2str(s),'s\n'])
    end
    if config.analysis
        m = floor(dt3/60);
        s = dt3-m*60;
        h = floor(m/60);
        m = m-h*60;
        fprintf(['分析总耗时',num2str(h),'h',num2str(m),'m',num2str(s),'s\n'])
    end
end

