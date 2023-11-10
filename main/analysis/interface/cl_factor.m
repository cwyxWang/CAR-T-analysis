function [cl,s,mask] = cl_factor(img_actin, s_mask, mask_index)

    s_real = s_mask==mask_index;
    mask = zeros(size(s_mask));
    if ~any(s_real(:))
        s = 0;
        cl = 0;
        return
    end
    
    s_edge = s_mask==-mask_index;

    CC = bwconncomp(s_real,26);
    s_label = labelmatrix(CC);
    cl = -inf;
    s_idx = 0;
    s = 0;
    for i = 1:CC.NumObjects
        s_seeds = s_edge;
        s_remain = s_label==i;
        dist_map = zeros(size(s_real));
        counter = 0;
        while any(s_seeds(:))
            temp = s_remain & imdilate(s_seeds,strel('sphere',1));
            s_seeds = temp;
            s_remain(temp) = 0;
            dist_map(temp) = counter;
            counter = counter+1;
        end
        max_dist = max(dist_map(:));
        if max_dist<1
            continue
        end
        y = double(img_actin(CC.PixelIdxList{i}));
        x = double(1-dist_map(CC.PixelIdxList{i})/max_dist);
        b = regress(y,[x,ones(size(x))]);
        if b(1)>cl
            cl = b(1);
            s_idx = i;
            s = size(CC.PixelIdxList{i},1);
        end
    end
    
    if s_idx==0
        cl=0;
    else
        mask(CC.PixelIdxList{s_idx}) = 1;
    end
    mask = mask>0;
end