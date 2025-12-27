% prueba minimización Least-Sq. con restricciones
close all
clear all
clc
load particle

%% 
sizec = size(C);
sized = size(d);

% variable de optimización x:
x = optimvar('x',sizec(2),'LowerBound',0);
% función objetivo
residual = C*x - d;
obj = sum(residual.^2);
% problema de optimización: 
nonneglsq = optimproblem('Objective',obj);
% optimoptions(nonneglsq) % imprime las opciones del problema
[sol,fval,exitflag,output] = solve(nonneglsq);
