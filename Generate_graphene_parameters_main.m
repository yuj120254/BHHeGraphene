%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% This file contains the code and examples to access the 5 different forms
% of the Graphene_paramters_function. 
% To run this file, you must have Graphene_lattice.m and
% Graphene_paramters_function.m in the same folder.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% The Graphene_paramters_function ran in this file can handle 5 different
% sets of inputs, each dealing with a different case of He on graphene, the
% cases are differentiated by the number of inputs, ranging from 3 to 7.
%
% The first 2 inputs does not change with the different cases and are:
% (strain type, array of strain indicies, .....)
% The strain type is either "isotropic" or "anisotropic"
% The array of strain indicies is an array of the indicies of the strain
% that you wish to generate graphene parameters for. To figure out which
% index corresponds to which strain, look at the params array below.
%
% The remaining inputs determine the case to calculate the parameters for
% Case 1: (mainly for testing)
% There are only 3 inputs: (strain type, strain indicies, "min")
% In this case, for each strain, the minimum of V(0,0,z) is determined, and
% the rest of the graphene parameters are calculated at that z.
%
% Case 2: (mainly for testing)
% There are 4 inputs: (strain type, strain indicies, "z", height(Angstroms))
% In this case, for each strain, the graphene parameters are calculated at
% the z specified.
%
% Case 3:
% There are 5 inputs: (strain type, strain indicies, "averaging", z, rho)
% In this case, the graphene parameters are calculated such that the
% potential is a weighted average in z using the weight rho. This case
% currently only supports calculations for a single value of strain.
%
% Case 4:
% There are 6 inputs: (strain type, strain indicies, "external_isotropic", V_min, V_sp, V_max)
% In this case, the graphene fourier components are calculated from the
% values of  V_min, V_sp, and V_max taken from some external potential, and
% the rest of the parameters are calculated from them. This case
% currently only supports calculations for a single value of strain.
%
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS DOES NOT CURRENTLY WORK
% Case 5: 
% There are 7 inputs: (strain type, strain indicies, "external_isotropic", V_min, V_sp_para, V_sp_perp, V_max)
% In this case, the graphene fourier components are calculated from the
% values of  V_min, V_sp, and V_max taken from some external potential, and
% the rest of the parameters are calculated from them. This case
% currently only supports calculations for a single value of strain.
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% THIS DOES NOT CURRENTLY WORK

% the paramters sigma and epsilon for all strain, using the latest values.
% The layout is:
% [[ strains .... ],
%  [ sigmas .... ],
%  [ epsilons .... ]]
anisotropic_params = [[0.00,0.05,0.10,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.34],
                      [2.643,2.669,2.700,2.738,2.749,2.759,2.769,2.780,2.791,2.804,2.817,2.831,2.846,2.858,2.878,2.896,2.917,2.938,3.065],
                      [16.961,16.758,16.530,16.256,16.170,16.086,16.001,15.911,15.825,15.720,15.618,15.510,15.400,15.303,15.156,15.025,14.882,14.733,13.963]];

isotropic_params = [[0.00,0.05,0.10,0.15,0.20,0.25,0.30],
                    [2.642,2.697,2.749,2.800,2.849,2.897,2.942],
                    [16.968,16.585,16.248,15.952,15.693,15.463,15.259]];
% Do not comment this out ^
                
             
    
% From here on uncomment the case which you need
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set the type of strain
% strain_type = "anisotropic";
% 
% % set the strains for which parameters are to be calclated
% is_array = [1,2,3,4,9,14];%
% 
% % display the strains for which parameters are to be calclated
% disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
% for s = is_array
%     strain = anisotropic_params(1, s)
% end
% 
% % Calls the functions
% [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%     Graphene_parameters_function(strain_type, is_array, "min");

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Case 1 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set the type of strain
% strain_type = "anisotropic";
% 
% % set the strains for which parameters are to be calclated
% is_array = [1,2,3,4,9,14];%
% 
% % set the height at which parameters are to be calclated (must be the same
% % lengths as is_array)
% z_array = [2.5,2.5,2.5,2.5,2.5,2.5];%
% 
% % display the strains for which parameters are to be calclated
% disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
% for s = is_array
%     strain = anisotropic_params(1, s)
% end
% 
% % Calls the functions
% [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%     Graphene_parameters_function(strain_type, is_array, "z", z_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Case 2 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set the type of strain
% strain_type = "isotropic";
% 
% % set the strains for which parameters are to be calclated
% is_array = [1,2,3,4];%
% 
% % load density in z for V0
% z_rho_V0 = dlmread("V0_rho0.txt");
% z_V0 = z_rho_V0(200:5:350,1);
% rho_V0 = z_rho_V0(200:5:350,2)/sum(z_rho_V0(200:5:350,2));% Important: always normalize rho here
% 
% z_array = [z_V0; z_V0; z_V0; z_V0];
% rho_array = [rho_V0; rho_V0; rho_V0; rho_V0];
% 
% % display the strains for which parameters are to be calclated
% disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
% for s = is_array
%     strain = anisotropic_params(1, s)
% end
% 
% % Calls the functions
% [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%     Graphene_parameters_function(strain_type, is_array, "average", z_array, rho_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Case 3 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% % set the type of strain
% strain_type = "isotropic";
% 
% % set the strains for which parameters are to be calclated
% is_array = [1];%
% 
% % load density in z for V0
% min_array = [-188.75805524953896];
% sp_array = [-160.58753714298558];
% max_array = [-154.69356294952843];
% 
% % display the strains for which parameters are to be calclated
% disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
% for s = is_array
%     strain = anisotropic_params(1, s)
% end
% 
% % Calls the functions
% [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%     Graphene_parameters_function(strain_type, is_array, "external_isotropic", min_array, sp_array, max_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Case 4 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%





%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% Case 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS DOES NOT CURRENTLY WORK

% % set the type of strain
% strain_type = "anisotropic";
% 
% % set the strains for which parameters are to be calclated
% is_array = [1];%
% % is_array = [1, 2, 19];%
% 
% % % load density in z for V0
% % min_array = [-188.75805524953896, -21.12602520878967];
% % sp_para_array = [-175.39950844354718, -7.98779019815916]; 
% % sp_perp_array = [-175.40070201328695, -7.987790198159161];
% % half_max_perp_array = [-171.61616608330277, -4.230340547886897];
% % max_para_array = [-171.61738906002037, -4.230340547886895];
% 
% % load density in z for V0
% min_array = [-188.75805524953896];
% sp_para_array = [-121.81671668775209]; 
% sp_perp_array = [-144.26849027314188];
% half_max_perp_array = [-143.55037721404412];
% max_para_array = [-153.75926337593884];
% 
% % % load density in z for V0
% % min_array = [-188.75805524953896, -184.19571474963163, -163.04332252894457];
% % sp_para_array = [-175.39950844354718, -168.29101289064602, -143.55037721404412]; 
% % sp_perp_array = [-175.40070201328695, -169.5638502855556, -153.75926337593884];
% % half_max_perp_array = [-171.61616608330277, -164.5929127221803, -145.2393372248427];
% % max_para_array = [-171.61738906002037, -163.61725278338128, 71.38735905486963];
% 
% % display the strains for which parameters are to be calclated
% disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
% for s = is_array
%     strain = anisotropic_params(1, s)
% end
% 
% % Calls the functions
% [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%     Graphene_parameters_function(strain_type, is_array, "external_anisotropic", ...
%     min_array, sp_para_array, sp_perp_array, half_max_perp_array, max_para_array);

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% End Case 5 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%THIS DOES NOT CURRENTLY WORK


