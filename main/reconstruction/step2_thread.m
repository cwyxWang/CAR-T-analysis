function timepoint = step2_thread(config)
% mip_num=校准体积数，mean噪声均值，std标准差，flipz偶数体积翻转，volume_num体积采样数
    if config.timepoint>=config.actual_test_num   
        timepoint = 0;
        return
    end
    
    config.h_info = get_data_info(config, 'h');
    config.k_info = get_data_info(config, 'k');
    
    timepoint = max(config.h_info.timepoint, config.k_info.timepoint);
    
    parfor i = config.actual_pass+1:timepoint
        
%         temp = floor((i-1)/config.volume_num);
%         volum_number = mod(temp,config.roi_num);
%         if volum_number~=18
%             continue
%         end

        tic;
        fprintf([num2str(i),'\t\t']);   % 打印时间点
        [save_rootdir, name_list] = arrange_data(i, config);
        
        step3_roi(config, i+config.timepoint, save_rootdir, name_list);
        
        dt = toc;
        fprintf([num2str(dt), '\n']);   % 打印时间间隔
    end

end
%%
function datainfo = get_data_info(config, type)

if type == 'h'
    data_path = config.h_path;
else
    data_path = config.k_path;
end
test_num = config.actual_test_num - config.timepoint;

datainfo = struct();
datainfo.data_path = data_path;
datainfo.type = type;
datainfo.data_name = {};
datainfo.stack_size_list = 0;
datainfo.width = 0;
datainfo.height = 0;
datainfo.bitdepth = 16;
datainfo.file_num = 0;
datainfo.timepoint = 0;

if isempty(data_path)
    return
end

if isfolder(data_path)
    temp = dir(fullfile(data_path,'*.tif'));   % *通配符，代表读取所有以.tif结尾的文件
    data_name = natsortfiles({temp.name});
    data_num = size(data_name, 2);  % data_num：有几张tif图
    stack_size_list = zeros(data_num,1);
    for i = 1:data_num
        temp = imfinfo(fullfile(data_path, data_name{i}));
        stack_size_list(i) = size(temp,1);   % 如有4张tif，则该列表保存了每张tif在z上的张数
        file_num = sum(stack_size_list);   % 打断，加速读取图像，在校准体积数为无穷时不起作用
        if file_num>test_num * config.slice_per_stack
            break
        end
    end
    info = temp(1);
    height   = info.Height;
    width    = info.Width;
    bitdepth = info.BitDepth;
else
    dcim = dcimg(data_path);
    data_name = {};
    stack_size_list = dcim.dc_sess_header.nfrms;
    width = dcim.dc_sess_header.xsize;
    height = dcim.dc_sess_header.ysize;
    bitdepth = 16;
end

file_num = sum(stack_size_list);
timepoint = ceil(file_num / config.slice_per_stack);   % floor向下取整
timepoint = min([timepoint,test_num]); 


datainfo.data_name = data_name;
datainfo.stack_size_list = stack_size_list;
datainfo.width = width;
datainfo.height = height;
datainfo.bitdepth = bitdepth;
datainfo.file_num = file_num;
datainfo.timepoint = timepoint;


end

%%
function [save_rootdir, name_list] = arrange_data(i, config)
%ARRANGE_DATA 此处显示有关此函数的摘要
%   此处显示详细说明
temp = floor((i-1)/config.volume_num);   % 计算出是哪一个ROI的一组图↓
volum_number = mod(temp,config.roi_num);   % mod取余，此部分是为了将（4）张tif里所有的slice进行重排，如第一个240张和第7个240张为一个细胞的图像
save_rootdir = fullfile(config.view_path,sprintf('\\MIP_ROI_%d\\',volum_number));
if ~exist(save_rootdir, 'dir')
    mkdir(save_rootdir)
end

a = (i-1) * config.slice_per_stack + 1;   % 取第i*60张图出来
b = i * config.slice_per_stack;
if config.flipz && mod(i,2) == 0   % 偶数翻转过程
    name_list = (b:-1:a);
else
    name_list = (a:b);
end

end





