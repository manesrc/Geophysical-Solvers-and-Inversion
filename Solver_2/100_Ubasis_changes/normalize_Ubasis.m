function norm_Ubasis = normalize_Ubasis(Ubasis)
% OUTPUT: same velocity basis, normalized by the velocity field with the highest norm
% INPUT: any velocity basis.
    norm_velo_max = 0;
    for i = 1:size(Ubasis,2)
        velo_i = Ubasis(:,i);
        norm_i = norm(velo_i);
        if norm_i > norm_velo_max
            norm_velo_max = norm_i;
        end 
    end 
    norm_Ubasis = Ubasis/norm_velo_max;
end 