function v = synapse_tracking(actin_mask, t_mask, s_mask)
if ~any(s_mask(:))
    v=0;
    return
end
t = size(t_mask,4);
t_center = zeros(t,3);
t_sum = zeros(1,3);
t_num = 0;
s_center = zeros(t,3);
s_sum = zeros(1,3);
s_num = 0;
for i = 1:t
    L = t_mask(:,:,:,i);
    stats = regionprops3(L>0);
    centers = stats.Centroid;
    vol = stats.Volume;
    if ~isempty(vol)
        [~,idx] = max(vol);
        t_center(i,:) = centers(idx,:);
        t_sum = t_sum + centers(idx,:);
        t_num = t_num + 1;
    end
end
t_mean = t_sum/t_num;
for i = 1:t
    L = s_mask(:,:,:,i);
    stats = regionprops3(L>0);
    centers = stats.Centroid;
    vol = stats.Volume;
    if ~isempty(vol)
        [~,idx] = max(vol);
        s_center(i,:) = centers(idx,:);
        s_sum = s_sum + centers(idx,:);
        s_num = s_num + 1;
    end
end
s_mean = s_sum/s_num;


[~, v_result] = tracker(actin_mask);
trace_num = size(v_result,1);
for i = 2:t
    if norm(t_center(i-1,:))>0 && norm(t_center(i,:))>0
        td = t_center(i,:) - t_center(i-1,:);
        if norm(td)<10
            for j = 1:trace_num
                if norm(squeeze(v_result(j,i,:)))>0
                    v_result(j,i,:) = squeeze(v_result(j,i,:))'-td;
                end
            end
        else    
            v_result(:,:,:) = 0.01;
            break
        end
    end
end


v_result = squeeze(sum(v_result,2));
if size(v_result,2)==1
    v_result = v_result';
end

vnorm = zeros(1,trace_num);
for i = 1:trace_num
    vnorm(i) = norm(v_result(i,:));
end
[~,idx] = sort(vnorm,'descend');
num = 4;
if size(idx,2)>=num
    s = zeros(1,3);
    for i = 1:num
        s = s + v_result(idx(i),:);
    end
    v_max = s/num;
else
    v_max = ones(1,3)*0.01;
end
% if isempty(v_result)
%     v_max = zeros(1,3);
% else
%     v_max = mean(v_result,1);
% end
if t_num>0 && s_num>0
    temp = t_mean - s_mean;
    temp = temp/norm(temp);
    v = temp*v_max';
else
    v = 0;
end

% v_all = zeros(1,t);
% for i = 2:t
%     if norm(s_center(i,:))>0 && norm(t_center(i,:))>0
%         cd = t_center(i,:) - s_center(i,:);
%         cd = cd/norm(cd);
%         s = zeros(1,trace_num);
%         for j = 1:trace_num
%             if norm(squeeze(v_result(j,i,:)))>0
%                 s(j) = cd*squeeze(v_result(j,i,:));
%             end
%         end
%         temp = s(s~=0);
%         if size(temp,2)>1
%             temp = sort(temp,'ComparisonMethod','abs');
%             v_all(i) = temp(end-1);
%         end
%     end
% end
% temp = v_all(v_all~=0);
% if isempty(temp)
%     v = 0;
% else
%     v = median(temp);
% end
end

function [trace_result, v_result] = tracker(actin_mask)
    timepoint = size(actin_mask,4);
    points = cell(1,timepoint);
    for i = 1:timepoint
        L = actin_mask(:,:,:,i);
        stats = regionprops3(L>0);
        centers = stats.Centroid;
        v = stats.Volume;
        data_point = cat(2,centers,v);
        points{i} = data_point;
    end

    min_volume = 0;
    max_volume = 500000;
    max_distance = 6;
    miss_tolerance = 0;
    min_length = 3;
    
    empty_trace = zeros(1,timepoint+3);
    empty_v = zeros(1,timepoint+3,3);
    trace_result = [];
    v_result = [];
    for i = 1:timepoint
        data_point = points{i};
        if isempty(trace_result)
            for j = 1:size(data_point,1)
                if data_point(j,4)>min_volume
                    new_trace = empty_trace;
                    new_trace(i) = j;
                    new_trace(i+1:i+3) = data_point(j,1:3);
                    trace_result = cat(1,trace_result,new_trace);
                    v_result = cat(1,v_result,empty_v);
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
                pend = detections(assignments(j,2), 1:3);
                pstart = trace_result(id, i:i+2);
                v_result(id, i, :) = pend-pstart;
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
                v_result = cat(1,v_result,empty_v);
            end
            
        end
    end
    
    i = 1;
    if ~isempty(trace_result)
        trace_result = trace_result(:,1:timepoint);
        v_result = v_result(:,1:timepoint,:);
        while i<=size(trace_result,1)
            if sum(trace_result(i,:)>0) < min_length
                trace_result(i,:) = [];
                v_result(i,:,:) = [];
            else
                i = i + 1;
            end
        end
    end
end
