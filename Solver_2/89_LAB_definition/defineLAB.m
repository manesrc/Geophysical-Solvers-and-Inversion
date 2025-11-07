function InfoLAB = defineLAB(curv_desc, varargin)
% defineLAB - Defines the X,Y coordinates of the LAB geometry.
%
% Syntax: InfoLAB = defineLAB(curv_desc, varargin)
%
% Inputs:
%   curv_desc (string): A case name specifying the geometry to use.
%       Examples: 'linear', 'sinusoidal', 'matrix_data', 'Afonso08'
%
%   varargin (cell array): Contains a structure with parameters required 
%       for the chosen 'curv_desc'.
%       - varargin{1}.model_width: Total width of the model [m]
%       - varargin{1}.model_height: Total depth of the model [m]
%
%       Case-specific inputs (passed inside varargin{1}):
%       - 'linear': 
%           - .mx: slope of the line
%           - .y0: y-intercept
%       - 'sinusoidal': (No additional inputs, uses hardcoded polynomial)
%       - 'matrix_data':
%           - .data_name: string, filename of the .mat file to load
%       - 'Afonso08': (No additional inputs, hardcoded)
%
% Outputs:
%   InfoLAB (structure): Contains the LAB coordinates
%       - .LABx: [1 x N] vector of X-coordinates [m]
%       - .LABy: [1 x N] vector of Y-coordinates [m]

    switch curv_desc
        case 'linear'
            % --- Inputs required:
            % varargin{1}.mx
            % varargin{1}.y0
            % varargin{1}.model_height
            % varargin{1}.model_width
            mx = varargin{1}.mx; 
            model_height = varargin{1}.model_height;
            model_width = varargin{1}.model_width; 
            y1 =  varargin{1}.y0; 
            LABy_f = @(x) y1 + mx.*x;
            InfoLAB.LABx = model_width * (0:0.1:1);    % [m]
            InfoLAB.LABy = LABy_f(InfoLAB.LABx); % [m] (dist. from surface)
        
        case 'sinusoidal'
            % --- Inputs required:
            % varargin{1}.model_height
            % varargin{1}.model_width (must be 660*1000)
            model_height = varargin{1}.model_height;
            model_width = varargin{1}.model_width;
            assert(model_width == 660*1000) % Check for specific model width
            % Solves a 3rd-order polynomial to fit points
            matrix_X = [0 0 0 1; 165^3 165^2 165 1; 
                495^3 495^2 495 1; 660^3 660^2 660 1];
            Y_points1 = [550; 600; 500; 550]-660;
            coef_curva = matrix_X\Y_points1;
            interphase = @(x) coef_curva(1)*x.^3 + coef_curva(2)*x.^2 + 1*coef_curva(3)*x + coef_curva(4);
            length1 = 30;
            InfoLAB.LABx = linspace(0,model_width,length1)';
            InfoLAB.LABy = model_height + 1000*interphase(InfoLAB.LABx/1000);
        
        case 'matrix_data'
            % --- Inputs required:
            % varargin{1}.data_name (e.g., 'my_lab_data.mat')
            % varargin{1}.model_height
            str_name = varargin{1}.data_name;
            a1 = load(str_name);
            % Assumes the .mat file contains 'datos_X' and 'datos_Y'
            InfoLAB.LABx = a1.datos_X*1000; % Convert km to m
            InfoLAB.LABy = varargin{1}.model_height + a1.datos_Y*1000; % Convert km to m
        
        case 'Afonso08'
            % --- Inputs required: None (hardcoded geometry)
            % This is a piecewise linear profile based on Afonso (2008)
            InfoLAB.LABx = [0 900 1010 1015 1050 1085 1915 1950 1985 1990 2100 3000]*1000; % [m]
            InfoLAB.LABy = [275 275 250 240 175 115 115 175 240 250 275 275]*1000; % [m]
    end 
end