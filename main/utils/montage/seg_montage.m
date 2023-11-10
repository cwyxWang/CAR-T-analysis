function seg_montage(view_path)

temp = dir(fullfile(view_path,'MIP_ROI_*'));
MIP_name = {temp.name};
mip_num = size(MIP_name, 2);

start = ceil(mip_num^0.5);
factor = mip_num / start;
while round(factor) ~= factor
    start = start+1;
    factor = mip_num / start;
end
row = factor;
colume = start;

bounder_width = 5;
center_line = 10;

mip_path_1 = fullfile(view_path,MIP_name{1},'ROIauto');
temp = dir(fullfile(mip_path_1,'*_predictions.h5'));
Img_name = {temp.name};
time_point = size(Img_name, 2);
h5_info = h5info(fullfile(mip_path_1,Img_name{1}), '/predictions');
size_data = h5_info.Dataspace.Size;

width = size_data(2);
height = size_data(1);
deepth = size_data(3);

xy_width = height*colume + bounder_width*(colume+1);
xy_height = (width+deepth)*row + bounder_width*(row*2+1)+center_line;

img_result = uint8(ones(xy_height, xy_width, time_point)*255);

for i = 1:mip_num
    
    colume_num = mod(i-1,colume);
    row_num = floor((i-1)/colume);
    
    w_start = colume_num*height + (colume_num+1)*bounder_width+1;
    w_end = w_start+height-1;
    
    h_start = row_num*width + (row_num+1)*bounder_width+1;
    h_end = h_start+width-1;
    
    z_start = row_num*deepth + row*width + (row+row_num+1)*bounder_width +center_line+ 1;
    z_end = z_start+deepth-1;
    
    mip_path_1 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'ROIauto');
    temp = dir(fullfile(mip_path_1,'*_predictions.h5'));
    Img_name = {temp.name};
    for j = 1:time_point
        try
            labels = h5read(fullfile(mip_path_1, Img_name{j}), '/predictions');
            % labels = labels(:,:,:,1);
            labels = 5-labels;
            labels(labels==5) = 0;
            img_result(h_start:h_end,w_start:w_end,j) = rot90(max(labels,[],3)*50);
            img_result(z_start:z_end,w_start:w_end,j) = flipud(permute(max(labels,[],2),[3,1,2])*50);
        catch
        end
    end
end
save_path = fullfile(view_path,'SegMontage');
if ~isfolder(save_path)
    mkdir(save_path);
end

for i = 1:time_point
    if i == 1
        imwrite(img_result(:,:,i),fullfile(save_path,'result.tif'));
%         imwrite(img_result(:,:,1,i),fullfile(save_path,'t.tif'));
%         imwrite(img_result(:,:,2,i),fullfile(save_path,'b.tif'));
%         imwrite(img_result(:,:,3,i),fullfile(save_path,'actin.tif'));
    else
        imwrite(img_result(:,:,i),fullfile(save_path,'result.tif'),'WriteMode','append');
%         imwrite(img_result(:,:,1,i),fullfile(save_path,'t.tif'),'WriteMode','append');
%         imwrite(img_result(:,:,2,i),fullfile(save_path,'b.tif'),'WriteMode','append');
%         imwrite(img_result(:,:,3,i),fullfile(save_path,'actin.tif'),'WriteMode','append');
    end
end

end

