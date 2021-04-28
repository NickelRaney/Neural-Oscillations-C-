function [poincare_map, point_time] = poincare_map(trajectory, plane, height)
% trajectory: (n, 3)
% plane: (2, 3)
orth_proj = orth(plane');
norm_vec = null(plane);
heights = trajectory * norm_vec;
poincare_map = [];
time = 0;
point_time = [];
for i =1:size(heights,1)
    if i>1
        if (heights(i)-height) * (heights(i-1)-height) <=0
            ratio = (heights(i-1) -height)/(heights(i-1) - heights(i));
            source = trajectory(i-1,:) * orth_proj;
            target = trajectory(i,:) * orth_proj;
            poincare_map = [poincare_map; source*ratio + target*(1-ratio)];
            if heights(i)-height <= 0
                point_time = [point_time; time];
            else
                point_time = [point_time; time];
                time = time + 1;
            end
        end
    end
end
end

