function write_test_data(labels, t_trace, save_path)
% labels = uint8(labels);
trace_num = size(t_trace,1);
fprintf('%d\n', trace_num);
t_name = fullfile(save_path,'t_test.tif');
b_name = fullfile(save_path,'b_test.tif');
d_name = fullfile(save_path,'d_test.tif');
surface_name = fullfile(save_path,'sur_test.tif');

for i = 1:size(t_trace,2)
    for k = 1:size(t_trace,3)
        i_label = labels(:,:,:,1,i,k);
        new_label = zeros(size(i_label));
        new_surface = new_label;
        for j = 1:trace_num
            if t_trace(j,i,k)>0
                new_label(-i_label==t_trace(j,i,k)) = j;
                new_surface(labels(:,:,:,4,i,k)==t_trace(j,i,k)) = j;
            end
        end
        labels(:,:,:,1,i,k) = new_label;
        labels(:,:,:,4,i,k) = new_surface;
    end
end
labels = uint8(labels);
for i = 1:size(labels,6)
    for j = 1:size(labels,5)
        img = max(labels(:,:,:,1,j,i),[],3);
        img2 = max(labels(:,:,:,2,j,i),[],3);
        img3 = max(labels(:,:,:,3,j,i),[],3);
        img4 = max(labels(:,:,:,4,j,i),[],3);
        if i == 1 && j==1
            imwrite(img, t_name);
            imwrite(img2, b_name);
            imwrite(img3, d_name);
            imwrite(img4, surface_name);
        else
            imwrite(img, t_name,'WriteMode','append');
            imwrite(img2, b_name,'WriteMode','append');
            imwrite(img3, d_name,'WriteMode','append');
            imwrite(img4, surface_name,'WriteMode','append');
        end
    end
end


end

