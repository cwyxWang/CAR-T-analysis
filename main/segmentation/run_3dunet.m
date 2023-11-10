function run_3dunet(data_path, ckpt, thread_ind, thread_num)
%RUN_3DUNET 此处显示有关此函数的摘要
%   此处显示详细说明
    gpu_ind = mod(thread_ind,thread_num);
    temp = dir(fullfile(data_path,'v_*'));
    View_name = {temp.name};
    view_num = size(View_name, 2);
    for i = 1:view_num
        view_path = fullfile(data_path,View_name{i});
        temp = dir(fullfile(view_path,'MIP_ROI_*'));
        Mip_name = natsortfiles({temp.name});
        mip_num = size(Mip_name, 2);
        for j = 1:mip_num
            if mod(j-1,thread_num) ~= thread_ind
                continue
            end
            mip_path = fullfile(view_path,Mip_name{j});
            temp = dir(fullfile(mip_path,'ROI*'));
            Roi_name = {temp.name};
            roi_num = size(Roi_name, 2);
            for k = 1:roi_num
                roi_path = fullfile(mip_path,Roi_name{k});
                temp = dir(fullfile(roi_path,'*.h5'));
                H5_name = {temp.name};
                h5_num = size(H5_name, 2);
                if h5_num>0
                    h5_info = h5info(fullfile(roi_path,H5_name{1}), '/raw');
                    data_size = h5_info.Dataspace.Size(1:3);
                    config_name = set_test_config(roi_path, ckpt, data_size);
                    system(sprintf('set CUDA_VISIBLE_DEVICES=%s && predict3dunet --config %s', num2str(gpu_ind), config_name));
                end
            end
        end
        seg_montage(view_path)
    end
    
    
    
end

