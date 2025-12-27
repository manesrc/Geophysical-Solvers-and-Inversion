function [final_velo, sol, fval2, exitflag, output, obj_val, optimName] = run_matlab_minim_vect_new_velo(alpha01, alpha02, nonlin, cond_dir, U_basis, InfoMesh, InfoProblem, InfoMaterial,Temp_est,pres_est)
    % this function computes difference between the temperature evaluated at the interface and the 
    % T_LAB isotherm to minimize it using the velocity basis. As the interface is frozen, the 
    % conductivity matrix should not change, therefore to avoid unnecessary computations it 
    % is calculated once.
    
    % set the temperature problem
    X_v = InfoMesh.X_v; 
    T_LAB = InfoProblem.T_LAB; 
    %T_ref = 0.005*T_LAB; 
    T_ref = T_LAB; 
    buffer_thickness = 80*1000; % 20 km for the gradual attenuation of convection velocities

    % Compute initial error
    crossed_elements = InfoMesh.list_cut; 
    % Temperature at Gamma_LAB
    figure_points = 0;
    % this computes many points at the interface 
    points2add = 3; 
    [MassMatrix,l_e] = compute_shapeFunctions_atLAB_morePoints(points2add, crossed_elements, figure_points, InfoMesh);
    num_cond_LAB = size(MassMatrix,1);
    cond_vect = T_LAB * ones(num_cond_LAB,1);
    % divide by the interface length in the element
    l_tot = sum(l_e); 
    
    % compute the attenuation factor gnl
    if nonlin == 1
        gnl = compute_g_nonlin(buffer_thickness,InfoMesh);
    else
        gnl = ones(size(MassMatrix,1),1); 
    end

    % objective function
    objective_function = @(alpha2)  (1/sqrt(l_tot)) * ( ( ( (MassMatrix * solveTemperatureProblem(gnl .* reshape(U_basis*alpha2,size(X_v)),Temp_est,pres_est,cond_dir,InfoMesh,InfoProblem,InfoMaterial) - cond_vect) / T_ref) ) .* sqrt(l_e)); 
    
    obj1 = objective_function(alpha01)' * objective_function(alpha01);
    obj2 = objective_function(alpha02)' * objective_function(alpha02); 

    if obj1 <= obj2
        alpha0 = alpha01;
    else
        alpha0 = alpha02;
    end 

    warning('off')

    scaledTolerance = 1e-5; 
    optimName = ['RealProblem_numUnk_',num2str(length(alpha0)),'_lin.txt'];
    scale_factor = 1;

    options1 = optimoptions('fsolve', ...
        'Algorithm','levenberg-marquardt', ...
        'Display', 'none', ...
        'OutputFcn', @(x,optimValues,state) funcValueStopFcn_vect(x, optimValues, state, scaledTolerance, optimName, scale_factor), ...  % Custom function to stop at desired Resnorm
        'TypicalX', alpha0, ...
        'FiniteDifferenceStepSize',1e-2,...
        'StepTolerance', 1e-8,...
        'MaxIterations',15,...
        'UseParallel',true,...
        'FunctionTolerance', 1e-6); 
    % ,...        'UseParallel',true
    % Solve the optimization problem with fminunc
    PROBLEM1.objective = @(alpha1) objective_function(alpha1);
    PROBLEM1.x0 = alpha0; 
    PROBLEM1.options = options1;
    PROBLEM1.solver = 'fsolve';
    % 
    % % Debugging: Check the initial state and dimensions
    % disp('Size of alpha0:');
    % disp(size(alpha0));
    % 
    % % Check the objective function value and size
    % try
    obj_val = objective_function(alpha0);
    %disp('squarenorm of the objective function value at initial alpha0:');
    %disp(norm(obj_val)^2);
    % catch ME
    %     disp('Error in objective function:');
    %     disp(ME.message);
    % end

    % Run optimization
    try
        [sol, fval, exitflag, output] = fsolve(PROBLEM1);
    catch ME
        disp('Error during optimization:');
        disp(ME.message);
    end
    
    % 
    % 
    % if sum((sol-alpha0) == 0) == length(alpha0)
    %     warning('the minimization did not perform any changes to the initial guess')
    % end 
    warning('on')
    fval2 = fval; 
    final_velo = gnl .* reshape(U_basis*sol,size(X_v));
end 


function stop = funcValueStopFcn_vect(x, optimValues, state, scaledTolerance, optName, scale_factor)
    stop = false;  % Continue optimization by default
    if (optimValues.fval)'*optimValues.fval < scaledTolerance
        stop = true;  % Stop if function value f(x) < scaledTolerance
    end
    
    % % Open the file in append mode
    % fid = fopen(optName, 'a');
    % 
    % % Log data based on the state of optimization
    % switch state
    %     case 'init'
    %         % Initialize the log file with a header
    %         fprintf(fid, 'Iter & FuncCount & ||f(x)||^2 & First-Ord Opt & lambda & Norm of step \\ \\ \n');
    %     case 'iter'
    %         % Log the iteration number, function value, and norm of step
    %         fprintf(fid, '%d & %d & %.2e & %.2e & %.2e & %.2e \\\\ \n', ...
    %                 optimValues.iteration, ...
    %                 optimValues.funccount, ...
    %                 ((optimValues.fval)'*optimValues.fval) / scale_factor, ...
    %                 optimValues.firstorderopt, ...
    %                 optimValues.lambda, ...
    %                 optimValues.stepsize);
    %     case 'done'
    %         fprintf(fid, 'Optimization completed.\n');
    % end
    % 
    % % Close the file
    % fclose(fid);
end
