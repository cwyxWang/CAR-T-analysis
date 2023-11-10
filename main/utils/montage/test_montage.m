function test_montage(view_path)

temp = dir(fullfile(view_path,'test*'));
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

info = imfinfo(fullfile(view_path,MIP_name{1},'t_test.tif'));
time_point = size(info, 1);
info = info(1);
width = info.Width;
height = info.Height;

xy_width = height*colume + bounder_width*(colume+1);
xy_height = width*row + bounder_width*(row+1);

img_t = uint8(ones(xy_height, xy_width, time_point)*255);
img_b = img_t;
img_s = img_t;

for i = 1:mip_num
    colume_num = mod(i-1,colume);
    row_num = floor((i-1)/colume);
    
    w_start = colume_num*height + (colume_num+1)*bounder_width+1;
    w_end = w_start+height-1;
    
    h_start = row_num*width + (row_num+1)*bounder_width+1;
    h_end = h_start+width-1;
    
    
    impath_1 = fullfile(view_path,sprintf('test%d',i-1),'t_test.tif');
    img_t(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point)*255;
    
    impath_1 = fullfile(view_path,sprintf('test%d',i-1),'b_test.tif');
    img_b(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point)*255;

    impath_1 = fullfile(view_path,sprintf('test%d',i-1),'sur_test.tif');
    img_s(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point)*128;

end
save_path = fullfile(view_path,'MontageTest');
if ~isfolder(save_path)
    mkdir(save_path);
end

for i = 1:time_point
    if i == 1
        imwrite(img_t(:,:,i),fullfile(save_path,'t.tif'));
        imwrite(img_b(:,:,i),fullfile(save_path,'b.tif'));
        imwrite(img_s(:,:,i),fullfile(save_path,'s.tif'));
    else
        imwrite(img_t(:,:,i),fullfile(save_path,'t.tif'),'WriteMode','append');
        imwrite(img_b(:,:,i),fullfile(save_path,'b.tif'),'WriteMode','append');
        imwrite(img_s(:,:,i),fullfile(save_path,'s.tif'),'WriteMode','append');
    end
end

end


function stack = read_mip_seq(path,h,w,t)
    stack = zeros(h,w,t);
    for i = 1:t
        try
            img = imread(path,i)';
            stack(:,:,i) = img(end:-1:1,:);
        catch
        end
    end
end
