function k_gp = compute_k_gp(one_material, relev_data,k1, k2,LS,Temp, pressure,N,Te,nel_x,InfoMaterial)
% compute_k_gp: computes the thermal conductivity at the gauss points
% depending on if the material is variable or not
% OUTPUT:
%   k_gp        : thermal conductivity at the gauss points
% INPUT:
%   one_material : flag to indicate if there is only one material in the
%                  element (1) or two materials (0)
%   relev_data  : relevant data to compute the mean thermal conductivity
%                  when there are two materials in the element
%   k1          : thermal conductivity of material 1
%   k2          : thermal conductivity of material 2
%   LS          : level set values at the gauss points
%   Temp        : temperature at the nodes of the element
%   pressure    : pressure at the nodes of the element
%   N           : shape functions at the gauss points
%   Te          : connectivity of the element
%   nel_x       : number of elements in the x direction
%   InfoMaterial: structure that contains the material properties
if InfoMaterial.k_var == 0
    k_gp = k1 * (LS >= 0) + k2 * (LS < 0);
    if one_material == 0
        delta_xi = ( relev_data(4) - relev_data(2) );
        delta_eta = (relev_data(5) - relev_data(3) );
        slope = delta_eta / delta_xi;
        mean_xi = (delta_xi)/2;
        eta_int = slope*mean_xi + relev_data(3);
        k_mean = 2 / ( (1-eta_int)/k1 + (1+eta_int)/k2 );
        k_gp = k_mean * ones(length(LS),1);            
    end 
else
    %error('revisar')
    T_sup = min(Temp);
    % material properties definitions
    k_0 = InfoMaterial.k_0;
    a_potencia = InfoMaterial.a;
    gamma = InfoMaterial.gamma;
    K0 = InfoMaterial.K0;
    derK0 = InfoMaterial.derK0;
    
    % to compute k_rad:
    k_rad1 = InfoMaterial.k_rad1;
    k_rad2 = InfoMaterial.k_rad2;
    k_rad3 = InfoMaterial.k_rad3;
    k_rad4 = InfoMaterial.k_rad4;
    
    % alpha coefficient:
    a1 = InfoMaterial.alphacoeff1;
    b1 = InfoMaterial.alphacoeff2;
    c1 = InfoMaterial.alphacoeff3; 
    n_dof = length(Temp);
    assert(length(Te) == 4)
    k_el = zeros(length(Te),1);
    for ii = 1:2
        kk  = [Te(1) Te(2)];
        dof_col_ii = (kk(ii):(nel_x+1):n_dof)';
        temp_in_col = Temp(dof_col_ii);
        press_in_col = pressure(dof_col_ii);
        alpha_in_col = a1 + b1 * temp_in_col + c1 * temp_in_col.^-2;
        int_in_col =zeros(size(dof_col_ii));
        for jj = length(int_in_col)-1:-1:1
            slope = (alpha_in_col(jj)-alpha_in_col(jj+1))/(temp_in_col(jj) - temp_in_col(jj+1));
            dT = temp_in_col(jj) - temp_in_col(jj+1); 
            int_in_col(jj) = int_in_col(jj+1) + slope*dT;
        end 
        k_rad_in_col = k_rad1 + k_rad2 * temp_in_col + k_rad3 * temp_in_col.^2 + k_rad4 * temp_in_col.^3;
        k_in_col = k_0 * ((T_sup./temp_in_col).^a_potencia) .* exp(-(4*gamma + 1/3 ) .* int_in_col ) .* (1 + derK0*(press_in_col)/(1000*K0)) + k_rad_in_col;
        k_el(2*ii-1:2*ii) = k_in_col(1:2); % change the position of k_el because now is [node1 node4 node2 node 3]
    end 
    k_el_order = [k_el(1); k_el(3); k_el(4); k_el(2)];
    k_gp = N*k_el_order; 
end 

end
 