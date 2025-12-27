function temperature = func_temp(pos,mx,Y0)

% generate a fake temperature
x1 = pos(:,1);
y1 = pos(:,2);

size_pos = size(pos,1);
temperature = zeros(size_pos,1);

%mx = -0.1;
%Y0 = 0.85;


T_LAB = 1505;
T_sup = 293;

for ii = 1:size_pos
    x_int = x1(ii);
    y_int = y1(ii);
    Y_LAB = Y0 + mx*x_int;
    if y_int <= Y_LAB           % below the LAB
        temp2save = T_LAB + (Y_LAB-y_int) * 0.5 * 660; 
    else
        temp2save = T_sup + ( (T_LAB-T_sup) / (1-Y_LAB)) * (1-y_int)  ;
    end 
    temperature(ii) = temp2save;
end 

end 