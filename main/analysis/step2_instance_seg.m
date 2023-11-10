function [raw_rearrange, labels_re, t_trace_re] = step2_instance_seg(roi_path, H5_name, timepoint, sample_num)
%%
    tic;
    fprintf('load data\t');
    temp = strsplit(H5_name{1}, '_predictions.h5');
    raw_name = [temp{1}, '.h5'];
    file_name = fullfile(roi_path, raw_name);
    h5_info = h5info(file_name, '/raw');
    size_raw = h5_info.Dataspace.Size;
    if size(size_raw, 2)==3
        size_raw = [size_raw, 1];
    end
    raw = uint16(zeros([size_raw, timepoint]));
    label = uint8(zeros([size_raw(1:3), timepoint]));
    parfor i = 1:timepoint
        temp = strsplit(H5_name{i}, '_predictions.h5');
        raw_name = [temp{1}, '.h5'];
        label_name = H5_name{i};
        label(:,:,:,i) = h5read(fullfile(roi_path, label_name), '/predictions');
        raw(:,:,:,:,i) = h5read(fullfile(roi_path, raw_name), '/raw');
    end
    toc;
%%
    tic;
    fprintf('instance segmention\t');
    labels = int8(zeros([size_raw(1:3),4, timepoint]));
    fliter = fspecial3('gaussian', 11, 3.5);
    assert(mod(size(fliter,1),2)==1);
    parfor i = 1:timepoint
        raw_temp = raw(:,:,:,:,i);
        if mod(i,sample_num)==0
            I = labels(:,:,:,:,i);
            I(:,:,:,3) = spot_detection(raw_temp(:,:,:,3), 450, fliter, 0.6);
            labels(:,:,:,:,i) = I;
            continue
        end
        I = label(:,:,:,i);
        t_mask = I==1;
        b_mask = I==2;
        actin = I==3;
        mem = I==4;
        
        D = imgaussfilt3(raw_temp(:,:,:,1),3);
        d_mask = imclose(D>398,strel('sphere',3));
        d_mask = bwareaopen(d_mask,20,26);
        
        
        t_instance = instance_seg(t_mask, actin);
        
        img488 = raw_temp(:,:,:,2);
%         img_down = imresize3(img488*8,0.5,'linear');
%         img488 = imresize3(img_down,size(img488),'linear');
%         J = imtophat(img488,strel('sphere',2));
%         temp = bwareaopen(J>=20,100,26);
        temp = spot_detection(img488, 400, fliter, 0.5);
        mask = int8(ones(size(temp)));
        mask(temp>0) = -1;
        t_instance = t_instance.*mask;
%         fprintf([num2str(sum(t_instance(:))),'\n'])
        b_instance = instance_seg(b_mask, mem);
        
        surface = int8(zeros(size(t_instance)));
        
        
        se_r = 10;
        se_dilate = strel('sphere',se_r);
        b_mask = b_instance~=0;
        not_b_dilate = ~imdilate(b_mask,se_dilate);
        % not_b_dilate = bwdist(b_mask)>se_r;
        b_round = imdilate(b_mask,strel('sphere',3));
        max_index = max(t_instance(:));
        for j = 1:max_index
            t_actin = abs(t_instance)==j;
            torch = b_round & t_actin;
            if ~any(torch(:))
                continue
            end
            t_dilate = imdilate(t_actin,se_dilate);
            % t_dilate = bwdist(t_actin)>se_r;
            background = ~t_dilate & not_b_dilate;
            CC = bwconncomp(background,26);
            background = labelmatrix(CC)==1;
            
            map_origin = b_mask*3 + t_actin*2 + background;
            [~,IDX] = bwdist(map_origin>0);
            map_exp = map_origin(IDX);
            surface_mask = map_exp==2 & imdilate(map_exp==3,strel('sphere',1));
            edge_mask = surface_mask & imdilate(map_exp==1,strel('sphere',se_r/2));
            surface(surface_mask) = j;
            surface(edge_mask) = -j;
        end
        
        I = labels(:,:,:,:,i);
        I(:,:,:,1) = t_instance;
        I(:,:,:,2) = b_instance;
        I(:,:,:,3) = d_mask;
        I(:,:,:,4) = surface;
        labels(:,:,:,:,i) = I;
    end
    toc;
%%
    tic;
    fprintf('tracking\t');
    t_instance = abs(squeeze(labels(:,:,:,1,:)));
    t_trace = tracker(t_instance);
    group_num = ceil(timepoint/sample_num);
    raw_rearrange = uint16(zeros([size_raw(1:3), 2, sample_num, group_num]));
    labels_re = int8(zeros([size_raw(1:3), 4, sample_num, group_num]));
    t_trace_re = uint8(zeros([size(t_trace,1), sample_num, group_num]));
    for i = 1:group_num
        for j = 1:sample_num
            idx = j+(i-1)*sample_num;
            if idx <= timepoint
                raw_rearrange(:,:,:,:,j,i) = raw(:,:,:,2:3,idx);
                labels_re(:,:,:,:,j,i) = labels(:,:,:,:,idx);
                try
                    t_trace_re(:,j,i) = t_trace(:,idx);
                catch
                end
            end
        end
    end
    toc;
end

%%
function mask = spot_detection(img, thr, fliter, r_thr)
    img = double(img);
    fliter_size = size(fliter, 1);
    r = floor(fliter_size/2);
    [height, width, deepth] = size(img);
    mask = int8(zeros([height width deepth]));
    for i = r+1:deepth-r
        for j = r+1:width-r
            for k = r+1:height-r
                if img(k,j,i)<thr
                    continue
                end
                img_patch = img(k-r:k+r, j-r:j+r, i-r:i+r);
                R = corrcoef(img_patch, fliter);
                if R(1,2)>r_thr
                    mask(k,j,i) = 1;
                end
            end
        end
    end
end

%%
function result = instance_seg(nucleus, membrane)
    result = int8(zeros(size(nucleus)));
    nucleus = bwareaopen(nucleus,3000,26);
    CC = bwconncomp(nucleus,26);
    nucleus_instance = int8(labelmatrix(CC));
    max_index = size(CC.PixelIdxList,2);
    if max_index == 0
        if sum(membrane(:))>3000
            CC = bwconncomp(membrane,26);
            max_pixel = 0;
            max_idx = 0;
            for j = 1:CC.NumObjects
                pixel_num = size(CC.PixelIdxList{j},1);
                if pixel_num>max_pixel
                    max_pixel = pixel_num;
                    max_idx = j;
                end
            end
            result(CC.PixelIdxList{max_idx}) = 1;
            temp = imfill(result>0,26,'holes');
            result = int8(temp);
        end
        return
    elseif max_index == 1
        membrane_instance = membrane;
    else
        [~,IDX] = bwdist(nucleus);
        membrane_instance = nucleus_instance(IDX);
        membrane_instance(~membrane) = 0;
    end
    
    for i = 1:max_index
        membrane_mask = membrane_instance==i;
        nucleus_mask = nucleus_instance==i;
        temp = membrane_mask & imdilate(nucleus_mask,strel('sphere',25));
        if sum(temp(:))<1000
            continue
        end
        CC = bwconncomp(membrane_mask,26);
        max_pixel = 0;
        max_idx = 0;
        for j = 1:CC.NumObjects
            pixel_num = size(CC.PixelIdxList{j},1);
            if pixel_num>max_pixel
                max_pixel = pixel_num;
                max_idx = j;
            end
        end
        nucleus_mask(CC.PixelIdxList{max_idx}) = 1;
        temp = imfill(nucleus_mask,26,'holes');
        result(temp) = i;
    end
    
    if sum(result(:))==0
        if sum(membrane(:))>3000
            CC = bwconncomp(membrane,26);
            max_pixel = 0;
            max_idx = 0;
            for j = 1:CC.NumObjects
                pixel_num = size(CC.PixelIdxList{j},1);
                if pixel_num>max_pixel
                    max_pixel = pixel_num;
                    max_idx = j;
                end
            end
            result(CC.PixelIdxList{max_idx}) = 1;
            temp = imfill(result>0,26,'holes');
            result = int8(temp);
        end
    end
end

%%
function trace_result = tracker(label)
    timepoint = size(label,4);
    points = cell(1,timepoint);
    for i = 1:timepoint
        L = label(:,:,:,i);
        stats = regionprops3(L);
        centers = stats.Centroid;
        v = stats.Volume;
        data_point = cat(2,centers,v);
        points{i} = data_point;
    end

    min_volume = 0;
    max_volume = 500000;
    max_distance = 200;
    miss_tolerance = 10;
    min_length = 50;
    
    empty_trace = zeros(1,timepoint+3);
    trace_result = [];
    for i = 1:timepoint
        data_point = points{i};
        if isempty(trace_result)
            for j = 1:size(data_point,1)
                if data_point(j,4)>min_volume
                    new_trace = empty_trace;
                    new_trace(i) = j;
                    new_trace(i+1:i+3) = data_point(j,1:3);
                    trace_result = cat(1,trace_result,new_trace);
                end
            end
        else
            trace_id = [];
            predictions = [];
            volume_id = [];
            detections = [];
            for j = 1:size(trace_result,1)
                if trace_result(j,i)>0||trace_result(j,i+1)>0||trace_result(j,i+2)>0
                    flag = false;
                    try
                        flag = sum(trace_result(j,i-miss_tolerance-1:i-1)) == 0;
                    catch
                    end
                    
                    if flag
                        trace_result(j,i:i+2) = 0;
                    else
                        trace_id = cat(1, trace_id, j);
                        predictions = cat(1, predictions, trace_result(j,i:i+2));
                    end
                end
            end
            for j = 1:size(data_point,1)
                if data_point(j,4)>min_volume && data_point(j,4)<max_volume
                    volume_id = cat(1, volume_id, j);
                    detections = cat(1, detections, data_point(j,1:3));
                end
            end
            if size(predictions,1) == 0
                assignments = [];
                unassignedTracks = [];
                unassignedDetections = (1:size(detections,1))';
            elseif size(detections,1) == 0
                assignments = [];
                unassignedTracks = (1:size(predictions,1))';
                unassignedDetections = [];
            else
                cost = zeros(size(predictions,1),size(detections,1));
                for j = 1:size(predictions, 1)
                    diff = detections - repmat(predictions(j,:),[size(detections,1),1]);
                    cost(j, :) = sqrt(sum(diff .^ 2,2));
                end
                [assignments,unassignedTracks,unassignedDetections] = assignDetectionsToTracks(cost,max_distance/2);
            end
            for j = 1:size(assignments, 1)
                id = trace_id(assignments(j,1));
                trace_result(id, i) = volume_id(assignments(j,2));
                trace_result(id, i+1:i+3) = detections(assignments(j,2), 1:3);
            end
            for j = 1:size(unassignedTracks, 1)
                id = trace_id(unassignedTracks(j));
                trace_result(id, i+1:i+3) = trace_result(id, i:i+2);
                trace_result(id, i) = 0;
            end
            for j = 1:size(unassignedDetections, 1)
                new_trace = empty_trace;
                new_trace(i) = volume_id(unassignedDetections(j));
                new_trace(i+1:i+3) = detections(unassignedDetections(j), 1:3);
                trace_result = cat(1,trace_result,new_trace);
            end
            
        end
    end
    
    i = 1;
    if ~isempty(trace_result)
        trace_result = trace_result(:,1:timepoint);
        while i<=size(trace_result,1)
            if sum(trace_result(i,:)>0) < min_length
                trace_result(i,:) = [];
            else
                i = i + 1;
            end
        end
    end
end
