function [coordinate] = rand_coordinate_generate(point, distance)
% Description:
% find a point that is a certain distance from another point
% x-axis y-axis 扩展的步长
step = distance;
x = point(1)  : point(1) + step;
y = point(2)  : point(2) + step;
[X,Y] = meshgrid(x,y);
distance_vector = sqrt((X-point(1)).^2 + (Y-point(2)).^2);
[row, col] = find(distance_vector == distance);
% 符合条件的坐标
coordiante_vector = [x(row)',y(col)'];
% 选取偏中间符合要求的坐标，避免只在x-axis 和 y-axis上
coordinate = coordiante_vector(size(coordiante_vector,1)-1,1:end);
end