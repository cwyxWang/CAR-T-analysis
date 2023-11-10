function count = dead_count(s_mask, b_mask, d_mask)

count = 0;
if ~any(s_mask(:))
    return
end
max_index = max(b_mask(:));
for i = 1:max_index
    b_instance = b_mask==i;
    touch = s_mask & imdilate(b_instance,strel('sphere',1));
    if ~any(touch(:))
        continue
    end
    dead = b_instance & imdilate(d_mask,strel('sphere',10));
    if any(dead(:))
        count = count + 1;
    end
end


end

