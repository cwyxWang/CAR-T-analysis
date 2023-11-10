function theta = tublin_theta(img_tublin, t_mask, tublin_mask, s_mask)

stats = regionprops3(s_mask);
s_centers = stats.Centroid;
s_num = size(s_centers,1);
if s_num==0
    theta = 180;
    return
end

intensity = 0;
index = 0;
CC = bwconncomp(t_mask & tublin_mask,26);
if CC.NumObjects==0
    theta = 180;
    return
end
for i = 1:CC.NumObjects
    mean_intensity = mean(img_tublin(CC.PixelIdxList{i}));
    if mean_intensity>intensity
        intensity = mean_intensity;
        index = i;
    end
end
tublin_mask = zeros(size(t_mask));
tublin_mask(CC.PixelIdxList{index}) = 1;

stats = regionprops3(t_mask);
centers = stats.Centroid;
volumes = stats.Volume;
[~,index] = max(volumes);
t_center = centers(index,:);

stats = regionprops3(tublin_mask);
centers = stats.Centroid;
tublin_center = centers(1,:);

theta_all = zeros(s_num, 1);
for i = 1:s_num
    v1 = tublin_center-t_center;
    v2 = s_centers(i,:)-t_center;
    cos_theta = v1*v2'/norm(v1)/norm(v2);
    theta_all(i) = acos(cos_theta)/pi*180;
end
theta = min(theta_all);

end


