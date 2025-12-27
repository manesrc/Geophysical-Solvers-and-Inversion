function name = generateNames2save(configuration_struct)
    % Generates a name for saving inversion results based on input parameters
    LAB_def  = configuration_struct.LAB_definition; 
    data_LAB = LAB_def.case;
    switch data_LAB
        case 'linear'    
            ini_depth = LAB_def.depth_ini;
            fin_depth = LAB_def.depth_ini + LAB_def.mx * configuration_struct.Long_X;
            meanDepth_ini = (ini_depth + fin_depth)/2; 
            if (meanDepth_ini < 620*1000) && (meanDepth_ini > 590*1000)
                str1 = '_sup_';
            elseif (meanDepth_ini < 590*1000) && (meanDepth_ini > 560*1000)
                str1 = '_int_';
            else
                str1 = '_inf_';
            end 
        otherwise 
            str1 = '';
    end     
    name = ['method',num2str(configuration_struct.forward_method),'_',configuration_struct.inversion_case,'_iniLAB',data_LAB,str1,'_Lprop_',...
            num2str(configuration_struct.len_prop),'_SigD_',num2str(configuration_struct.sigma_d),'_nx_',num2str(configuration_struct.nel_x),'_ny_',num2str(configuration_struct.nel_y),...
            '_nst_',num2str(configuration_struct.nstep),'.mat'];
end 