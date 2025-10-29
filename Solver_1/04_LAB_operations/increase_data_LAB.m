function [datos_x,datos_y] = increase_data_LAB(InfoLAB,numpoints2add)
% this function will increase the number of data in the vectors 
%% OUTPUT 
% datos_x: new vector with linear interpolation to find new points [length(dataLAB_x)*(data_increase-2)]
% datos_y: new vector with linear interpolation to find new points [length(dataLAB_y)*(data_increase-2)]
%% INPUT
% InfoLAB: contains the data of the points in the LAB 
% data_increase: increasing amount of points minus two, the beginning and
% end of the original line
%% code
dataLABx = InfoLAB.LABx;
dataLABy = InfoLAB.LABy;
orig_length = length(dataLABx);

datos_x = zeros(orig_length*(numpoints2add-2),1);
datos_y = zeros(orig_length*(numpoints2add-2),1);

fin_index = 1;
for ii = 1:orig_length-1
    datos_x_ii = [dataLABx(ii) dataLABx(ii+1)];
    datos_y_ii = [dataLABy(ii) dataLABy(ii+1)];
    points_added_x = linspace(datos_x_ii(1), datos_x_ii(2),numpoints2add);
    points_added_y = interp1(datos_x_ii,datos_y_ii,points_added_x);
    ini_index = fin_index;
    fin_index = ini_index + (numpoints2add-1);
    datos_x(ini_index:fin_index) = points_added_x';
    datos_y(ini_index:fin_index) = points_added_y';
end 

end 