function [LABx,LABy] = defineLAB(informationLAB)
% this function defines the LAB disposition 
% OUTPUTS: 
% LABx: points belonging to the interface in X direction [from 0 to L_ref]
% LABy: points belonging to the interface in Y direction [in km in depth from the surface]
% INPUTS: 
% informationLAB: structure
%               % three cases available 1. Linear, 2. Sinusoidal, 3. data
%                if case linear needs point1 and point2 in y direction to draw a line

%% 
caso_aplicado = informationLAB.disposition;
datosLAB = open('datos_isotherm.mat');

switch(caso_aplicado)
    case 0 % linear distribution
        % toy case [horizontal/inclined LAB]
        length1 = length(datosLAB.datos_X)+1;
        LABx = linspace(0,660,length1)';
        LABy = -660*linspace(informationLAB.y_ini,informationLAB.y_fin,length1)';
    case 1 % sinusoidal
        % % resolviendo coef = X\Y da indeterminado y después arrastra mucho error
        matrix_X = [0 0 0 1; 165^3 165^2 165 1; 
            495^3 495^2 495 1; 660^3 660^2 660 1];
        Y_points1 = [550; 600; 500; 550]-660;
        coef_curva = matrix_X\Y_points1;
        interphase = @(x) coef_curva(1)*x.^3 + coef_curva(2)*x.^2 + 1*coef_curva(3)*x + coef_curva(4);
        length1 = length(datosLAB.datos_X)+1;
        LABx = linspace(0,660,length1)';
        LABy = interphase(LABx);
    case 2 % data
        LABx = flip(datosLAB.datos_X); 
        LABy = flip(datosLAB.datos_Y);
end 

