function make_montage(view_path)

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

mip_path_1 = fullfile(view_path,MIP_name{1},'MIPxy_405');
temp = dir(fullfile(mip_path_1,'*.tif'));
Img_name = {temp.name};
time_point = size(Img_name, 2);
info = imfinfo(fullfile(mip_path_1,Img_name{1}));
width = info.Width;
height = info.Height;

info = imfinfo(fullfile(view_path,MIP_name{1},'MIPxz_405',Img_name{1}));
deepth = info.Width;

xy_width = height*colume + bounder_width*(colume+1);
xy_height = (width+deepth)*row + bounder_width*(row*2+1)+center_line;

img_405 = uint16(zeros(xy_height, xy_width, time_point));
img_488 = img_405;
img_561 = img_405;

for i = 1:mip_num
    colume_num = mod(i-1,colume);
    row_num = floor((i-1)/colume);
    
    w_start = colume_num*height + (colume_num+1)*bounder_width+1;
    w_end = w_start+height-1;
    
    h_start = row_num*width + (row_num+1)*bounder_width+1;
    h_end = h_start+width-1;
    
    z_start = row_num*deepth + row*width + (row+row_num+1)*bounder_width +center_line+ 1;
    z_end = z_start+deepth-1;
    
    impath_1 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxy_405');
    impath_2 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxz_405');
    if isfolder(impath_1)
        img_405(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point);
        img_405(z_start:z_end,w_start:w_end,:) = read_mip_seq(impath_2, deepth, height, time_point);
        
    end
    impath_1 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxy_488');
    impath_2 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxz_488');
    if isfolder(impath_1)
        img_488(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point);
        img_488(z_start:z_end,w_start:w_end,:) = read_mip_seq(impath_2, deepth, height, time_point);
    end
    impath_1 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxy_561');
    impath_2 = fullfile(view_path,sprintf('MIP_ROI_%d',i-1),'MIPxz_561');
    if isfolder(impath_1)
        img_561(h_start:h_end,w_start:w_end,:) = read_mip_seq(impath_1, width, height, time_point);
        img_561(z_start:z_end,w_start:w_end,:) = read_mip_seq(impath_2, deepth, height, time_point);
    end
end
save_path = fullfile(view_path,'Montage');
if ~isfolder(save_path)
    mkdir(save_path);
end

img_405(img_405==0) = max(img_405(:));
img_488(img_488==0) = max(img_488(:));
img_561(img_561==0) = max(img_561(:));
for i = 1:time_point
    if i == 1
        imwrite(img_405(:,:,i),fullfile(save_path,'405.tif'));
        imwrite(img_488(:,:,i),fullfile(save_path,'488.tif'));
        imwrite(img_561(:,:,i),fullfile(save_path,'561.tif'));
    else
        imwrite(img_405(:,:,i),fullfile(save_path,'405.tif'),'WriteMode','append');
        imwrite(img_488(:,:,i),fullfile(save_path,'488.tif'),'WriteMode','append');
        imwrite(img_561(:,:,i),fullfile(save_path,'561.tif'),'WriteMode','append');
    end
end

end


function stack = read_mip_seq(path,h,w,t)
    stack = zeros(h,w,t);
    temp = dir(fullfile(path,'*.tif'));
    Img_name = {temp.name};
    img_num = size(Img_name, 2);
    for i = 1:img_num
        img = imread(fullfile(path,Img_name{i}))';
        stack(:,:,i) = img(end:-1:1,:);
    end
end
