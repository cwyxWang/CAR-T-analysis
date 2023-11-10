function step3_analysis(raw, labels, t_trace, save_path, ind)
timepoint = size(t_trace,3);
sample_num = size(t_trace,2);
v_all = zeros(1, timepoint);
s_all = zeros(1, timepoint);
theta_all = zeros(1, timepoint);
d_all = zeros(1, timepoint);
cl_all = zeros(1, timepoint);

parfor i = 1:timepoint
%     fprintf([num2str(i),'\n']);
    img = raw(:,:,:,:,:,i);
    I = labels(:,:,:,:,:,i);
    tt = t_trace(1,:,i);
    
    img_tublin = img(:,:,:,2,sample_num);
    tublin_mask = I(:,:,:,3,sample_num);
    t_mask = zeros(size(tublin_mask))>0;
    s_mask = t_mask;
    for j = sample_num-1:-1:1
        if tt(j)>0
            mask_index = tt(j);
            t_mask = abs(I(:,:,:,1,j))==mask_index;
            s_mask = I(:,:,:,4,j)==mask_index;
            break
        end
    end
    theta_all(i) = tublin_theta(img_tublin, t_mask, tublin_mask, s_mask);
%     
%     t_mask = squeeze(I(:,:,:,1,:));
%     b_mask = squeeze(I(:,:,:,2,:));
%     d_mask = squeeze(I(:,:,:,3,:));
%     s_mask = squeeze(I(:,:,:,4,:));
    xyz = size(I);
    d_temp = zeros(1,sample_num-1);
    cl_temp = d_temp;
    s_temp = d_temp;
    s_tracked = zeros([xyz(1:3),sample_num-1]);
    t_sport = s_tracked;
    t_mask = s_tracked;
    for j = 1:sample_num-1
        mask_index = int8(tt(j));
        smaskj = I(:,:,:,4,j);
        bmaskj = I(:,:,:,2,j);
        dmaskj = I(:,:,:,3,j);
        if mask_index>0
            [cl, s, mask] = cl_factor(img(:,:,:,1,j), smaskj, mask_index);
            cl_temp(j) = cl;
            s_temp(j) = s;
            s_tracked(:,:,:,j) = mask;
            temp = I(:,:,:,1,j);
            t_sport(:,:,:,j) = temp==-mask_index;
            t_mask(:,:,:,j) = abs(temp)==mask_index;
            d_temp(j) = dead_count(smaskj==mask_index, bmaskj, dmaskj>0);
        end
    end
    d_all(i) = max(d_temp);
    if any(s_temp>0)
        cl_all(i) = mean(cl_temp(s_temp>0));
        s_all(i) = mean(s_temp(s_temp>0));
    else
        cl_all(i) = 0;
        s_all(i) = 0;
    end
    
    v_all(i) = synapse_tracking(t_sport>0, t_mask>0, s_tracked>0);
    
end

ind = num2str(ind);
% fprintf([ind,'\n']);
file_name  = fullfile(save_path,'actin_result.csv');
write_csv(file_name,v_all,ind,0);
file_name  = fullfile(save_path,'surface_result.csv');
write_csv(file_name,s_all,ind,0);
file_name  = fullfile(save_path,'claer_result.csv');
write_csv(file_name,cl_all,ind,0);
file_name  = fullfile(save_path,'kill_result.csv');
write_csv(file_name,d_all,ind,-1);
file_name  = fullfile(save_path,'tublin_result.csv');
write_csv(file_name,theta_all,ind,180);
end
%%
function write_csv(file_name,array,ind,void)
    fid = fopen(file_name, 'a+', 'n', 'utf8');
    len = size(array,2);
    str = ind;
    for i = 1:len
        if array(i)~=void
            str = [str,',',num2str(array(i))];
        else
            str = [str,','];
        end
    end
    str = [str,'\n'];
    fprintf(fid, str);
    fclose(fid);
end

