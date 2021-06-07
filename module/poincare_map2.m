function [poincare_map] = poincare_map2(trajectory, plane, height)
% trajectory: (n, 3)
% plane: (2, 3)
orth_proj = orth(plane');
norm_vec = null(plane);
heights = trajectory * norm_vec;
poincare_map = [];
time = 0;
point_time = [];

j=-100;

for i =21:size(heights,1)-20
    if (heights(i)-height) * (heights(i-1)-height) <=0 && (mean(heights(i-20:i-1))-height)<0 && (mean(heights(i:i+20))-height)>500 && i-j>200
        ratio = (heights(i-1) -height)/(heights(i-1) - heights(i));
        source = trajectory(i-1,:) * orth_proj;
        target = trajectory(i,:) * orth_proj;
        poincare_map = [poincare_map; source*ratio + target*(1-ratio)];
        j=i;
    end
end
end

