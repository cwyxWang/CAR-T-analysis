function step3_roi(config, index, save_rootdir, name_list)
%MIP_VOLUM 此处显示有关此函数的摘要
%   此处显示详细说明 ROIfull
    config.save_rootdir = save_rootdir;
    config.name_list = name_list;
    
    result = step4_reconstruct(config);
    
    roi_dir = fullfile(config.save_rootdir, 'ROIauto');
    if ~isfolder(roi_dir)
        mkdir(roi_dir)
    end
    
    if config.view_center == 0
        roi = result;
    elseif config.view_center > 0
        roi_center = config.view_center;
        r = floor(size(result, 1)/2);
        roi = result(:, roi_center-r:roi_center-r+size(result, 1)-1, :, :);
    else
        
        roi_text = fullfile(roi_dir, sprintf('roi_binning%d.txt', config.binning));
        if ~isfile(roi_text)
            mip_488 = max(result(:,:,:,2),[],3);
            roi_center = get_w_center(mip_488, config.binning);
            fid = fopen(roi_text,'a');
            fprintf(fid, '%d\n',roi_center);
            fclose(fid);
        end
        fid = fopen(roi_text);
        temp = textscan(fid,'%d');
        roi_center = temp{1}(1);
        fclose(fid);

        r = floor(size(result, 1)/2);
        roi = result(:, roi_center-r:roi_center-r+size(result, 1)-1, :, :);
    end
    
    write_roi(config, index, roi, roi_dir);   
end


%%
function center = get_w_center(img,binning)
    r = 280/binning;
    wr = 160/binning;
    img(:,end-wr:end) = 0;
    img1 = imresize(img,0.5);
    img2 = imresize(img1,0.5);
    img3 = imresize(img2,0.5);
    img4 = imresize(img3,0.5);
    center4 = serch_step(img4, ceil(r/16), 0, 0);
    center3 = serch_step(img3, ceil(r/8), center4, 1.5);
    center2 = serch_step(img2, ceil(r/4), center3, 1.3);
    center1 = serch_step(img1, ceil(r/2), center2, 1.2);
    center = serch_step(img, r, center1, 1.1);
end

function center = serch_step(img, r, center, serch_range)
    if center>0
        width = size(img, 2);
        lower_limit = max(1,floor(center*2-r*serch_range));
        higher_limit = min(ceil(center*2+r*serch_range),width);
        img = img(:,lower_limit:higher_limit);
    else
        lower_limit = 1;
    end
    se = strel('disk',r);
    height = size(img, 1);
    cut_off = ceil((2*r+1-height)/2);
    kernel = se.Neighborhood(cut_off:cut_off+height-1,:);
    score = conv2(img,kernel,'valid');
    [~,center] = max(score);
    center = center+r+lower_limit-1;
end

%%
function write_roi(config, index, roi, roi_dir)
%WRITE_STACK 此处显示有关此函数的摘要
%   此处显示详细说明
channel_num = size(roi,4);
input_channel = 3;
channel_name = cell(3,1);
channel_name{1} = '405';
channel_name{2} = '488';
channel_name{3} = '561';

for i = 1:channel_num
    stack = roi(:,:,:,i);
    mip_xz = squeeze(max(stack,[],2));
    mip_xy = squeeze(max(stack,[],3));
    
    if ~config.rotate
        mip_xy = imwarp(mip_xy, affine2d([config.factor*sin(config.theta) 0 0; 0 1 0; 0 0 1]));
    end
    
    mipxy_name = fullfile(config.save_rootdir,['MIPxy_',channel_name{i}]);
    if ~exist(mipxy_name, 'dir')
        mkdir(mipxy_name);
    end
    mipxz_name = fullfile(config.save_rootdir,['MIPxz_',channel_name{i}]);
    if ~exist(mipxz_name, 'dir')
        mkdir(mipxz_name);
    end
    
    imwrite(mip_xy, fullfile(mipxy_name, sprintf('%05d.tif', index)));
    imwrite(mip_xz, fullfile(mipxz_name, sprintf('%05d.tif', index)));
end

if config.h5
    file_name = fullfile(roi_dir, sprintf('%05d.h5', index));
    if exist(file_name, 'file')
        delete(file_name);
    end
    
    if channel_num>input_channel
        raw = roi(:,:,:,1:input_channel);
        raw_2 = roi(:,:,:,input_channel+1:end);
        h5create(file_name, "/raw",size(raw),'Datatype','uint16');
        h5write(file_name,"/raw",raw);
        h5create(file_name, "/raw_2",size(raw_2),'Datatype','uint16');
        h5write(file_name,"/raw_2",raw_2);
    elseif channel_num == input_channel
        h5create(file_name, "/raw",size(roi),'Datatype','uint16');
        h5write(file_name,"/raw",roi);
    else
        h5create(file_name, "/raw_2",size(roi),'Datatype','uint16');
        h5write(file_name,"/raw_2",roi);
    end
    
else
    for i = 1:channel_num
        new_folder = fullfile(roi_dir, channel_name{i});
        if ~exist(new_folder, 'dir')
            mkdir(new_folder)
        end
        file_name = fullfile(new_folder, sprintf('%05d.tif', index));
        deepth = size(roi,3);
        for j = 1:deepth
            if j == 1
                imwrite(roi(:,:,j,i), file_name);
            else
                imwrite(roi(:,:,j,i), file_name,'WriteMode','append');
            end
        end
    end
end

end