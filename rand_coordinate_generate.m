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
coordiante_vecotr = [x(row)',y(col)'];
coordinate = coordiante_vecotr(size(coordiante_vecotr,1)-1,1:end);
end