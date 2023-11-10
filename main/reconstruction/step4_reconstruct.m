function stack_result = step4_reconstruct(config)
%RECONSTRUCTION_STEP 此处显示有关此函数的摘要
%   此处显示详细说明
    shift_matrix = get_shift(config.factor, config.theta, config.new_system);
    rotate_matrix = get_rotate(config.factor, config.theta, config.new_system);
    % affine_matrix = get_affine(config.factor, config.theta, config.new_system, false);

    h_stack = read_stack(config.h_info, config.name_list, config.new_system, config.wview);
    k_stack = read_stack(config.k_info, config.name_list, config.new_system, false);
    k_stack = k_stack(end:-1:1,:,:);
    if ~isempty(h_stack) && ~isempty(k_stack) && (size(h_stack,2) ~= size(k_stack,2) || size(h_stack,1) ~= size(k_stack,1))
        warning('图像尺寸不一致，将进行裁剪')
        h_min = min(size(h_stack,1),size(k_stack,1));
        w_min = min(size(h_stack,2),size(k_stack,2));
        stack = cat(4, h_stack(1:h_min,1:w_min,:,:), k_stack(1:h_min,1:w_min,:));
    else
        stack = cat(4, h_stack, k_stack);
    end
    if config.binning>1
        stack = imresize(stack*config.binning^2,1/config.binning,'bilinear');
    end
    
    stack_result = [];
    
    for j = 1:size(stack,4)
        
        img_channel = stack(:,:,:,j);
        result = imwarp(img_channel,affine3d(shift_matrix), 'linear');         % 'cubic','nearest'
        if config.rotate
            result = imwarp(result,affine3d(rotate_matrix), 'linear');
            depth = floor(size(result, 3)/2);
            depth_goal = floor(size(stack,1)/2*sin(config.theta));
            result = result(:,:,depth-depth_goal+1:depth+depth_goal);
        end
        
        if config.mean>0
            noise = uint16(randn(size(result))*config.std + config.mean);
            flag_ = result<config.mean-3*config.std;
            result(flag_) = noise(flag_);
        end
        
        result = uint16(result);
        stack_result = cat(4,stack_result,result);
    end
end

%%
function shift_matrix = get_shift(factor, theta, new_system)
    if new_system
        a = factor * cos(theta);
    else
        a = -factor * cos(theta);
    end
    shift_matrix = [1           0   0            0
                    0           1   0            0
                    0           a   1            0
                    0           0   0            1];
end

function correction_matrix = get_rotate(factor, theta, new_system)
    c = cos(theta);
    s = sin(theta);
    if new_system
        correction_matrix = [0               1  0            0
                             c               0  -s           0
                             factor*s*s      0  factor*s*c   0
                             0               0  0            1];
    else
        correction_matrix = [0               1  0            0
                             -c              0  s            0
                             factor*s*s      0  factor*s*c   0
                             0               0  0            1];
    end
end

% function affine_matrix = get_affine(factor, theta, new_system, iso)
%     if iso
%         a = 1;
%     else
%         a = 1/sin(theta)/factor;
%     end
%     if new_system
%         affine_matrix = [0               1  0            0
%                          cos(theta)*a    0  -sin(theta)  0
%                          factor*a        0  0            0
%                          0               0  0            1];  
%     else
%         affine_matrix = [0               1  0            0
%                          -cos(theta)*a   0  sin(theta)   0
%                          factor*a        0  0            0
%                          0               0  0            1];
%     end
% end

%%
function stack = read_stack(iminfo, name_list, new_system, wview)
n_slice  = size(name_list, 2);   % name_list是一个1*60的矩阵，读取第二个维度60给n_slice
if isempty(iminfo.data_path)
    stack = [];
    return
end
if iminfo.type == 'h'
    name_list = name_list-1;
    name_list(name_list<1) = 1;
end
if isfolder(iminfo.data_path)
    name_code = cell(1,n_slice);   % cell元胞数组，存字符串的数组，索引用{}大括号
    slice_code = (1:n_slice)*0;
    img_num = size(iminfo.stack_size_list);
    
    for i = 1:n_slice
        for j = 1:img_num
            if name_list(i)<=sum(iminfo.stack_size_list(1:j))
                name_code{i} = iminfo.data_name{j};
                slice_code(i) = name_list(i)-sum(iminfo.stack_size_list(1:j-1));
                break
            end
        end
    end
    
    if iminfo.type == 'k'   % K开头的没有Wview，单通道，因此第四个维度为1
        if new_system
            stack = zeros(iminfo.height, iminfo.width, n_slice, 1);
        else
            stack = zeros(iminfo.width, iminfo.height, n_slice, 1);
        end
    elseif iminfo.bitdepth == 48 % 48位代表RGB
        stack = zeros(iminfo.height, iminfo.width, n_slice, 2);
    else   % 如果H开头有Wview，双通道，则第四个维度为2
        if wview
            h = floor(iminfo.height/2);   % 上下拆一半
            stack = zeros(h, iminfo.width, n_slice, 2);
        else
            stack = zeros(iminfo.height, iminfo.width, n_slice, 1);
        end
    end
    
    for i = 1:n_slice
        if iminfo.type == 'k'
            try
                if new_system
                    stack(end:-1:1,:,i,:) = imread(fullfile(iminfo.data_path, name_code{i}), slice_code(i));
                    % stack(:,:,i,:) = imread(fullfile(iminfo.data_path, name_code{i}), slice_code(i));
                else
                    stack(:,:,i,:) = imread(fullfile(iminfo.data_path, name_code{i}), slice_code(i))';
                end
            catch
                warning('索引超出数组元素的数目.');
            end
        elseif iminfo.bitdepth == 48
            try
                img_up = imread(fullfile(iminfo.data_path, name_code{i}), 1);
                img_down = imread(fullfile(iminfo.data_path, name_code{i}), 2);
                stack(:,:,i,:) = cat(4, img_up, img_down);
            catch
                warning('索引超出数组元素的数目.');
            end
        else
            if wview
                try
                    temp = imread(fullfile(iminfo.data_path, name_code{i}), slice_code(i));
                    img_up = temp(1:h,:);
                    img_down = temp(iminfo.height-h+1:iminfo.height,:);
                    stack(:,:,i,:) = cat(4, img_up, img_down);   % cat：在第四个维度上把上下拼接到一起
                catch
                    warning('索引超出数组元素的数目.');
                end
            else
                try
                    stack(:,:,i,:) = imread(fullfile(iminfo.data_path, name_code{i}), slice_code(i));
                catch
                    warning('索引超出数组元素的数目.');
                end
            end
        end
    end
else
    stack = zeros(iminfo.height/2, iminfo.width, n_slice, 2);
    for i = 1:n_slice
        try
            dcim = dcimg(iminfo.data_path, name_list(i));
            img_up = dcim.data(1:2:end,:,:);
            img_down = dcim.data(2:2:end,:,:);
            img_up(1,:) = 0;
            img_down(1,:) = 0;
            stack(:,:,i,:) = cat(4, img_up, img_down);
        catch
            warning('索引超出数组元素的数目.');
        end
    end
end
end
