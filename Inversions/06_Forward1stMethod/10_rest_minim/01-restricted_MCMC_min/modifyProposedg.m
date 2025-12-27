function new_g = modifyProposedg(sigma,old_g)

nDOF = length(old_g);

nmodif = 1;%ceil(nDOF/3);

params = 1:nDOF; % all the DOF
id_change_param =  randperm(nDOF,nmodif)  ; % Determine one to modify
p = params(id_change_param);

id_2change = zeros(nDOF,1);
id_2change(p) = 1;

g_change = id_2change.*(sigma*(1/3)*randn(nDOF,1));
new_g = old_g + g_change;

end
