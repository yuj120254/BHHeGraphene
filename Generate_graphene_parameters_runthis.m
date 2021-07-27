% the paramters sigma and epsilon for all strain, using the latest values
% layout is:
% [[ strains .... ],
%  [ sigmas .... ],
%  [ epsilons .... ]]
anisotropic_params = [[0.00,0.05,0.10,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.34],
                      [2.643,2.669,2.700,2.738,2.749,2.759,2.769,2.780,2.791,2.804,2.817,2.831,2.846,2.858,2.878,2.896,2.917,2.938,3.065],
                      [16.961,16.758,16.530,16.256,16.170,16.086,16.001,15.911,15.825,15.720,15.618,15.510,15.400,15.303,15.156,15.025,14.882,14.733,13.963]];

isotropic_params = [[0.00,0.05,0.10,0.15,0.20,0.25,0.30],
                    [2.642,2.697,2.749,2.800,2.849,2.897,2.942],
                    [16.968,16.585,16.248,15.952,15.693,15.463,15.259]];
    
% indicates the type of strain 
strain_type = "anisotropic";

% load density in z for V0
z_rho_V0 = dlmread("V0_rho0.txt");
z_V0 = z_rho_V0(200:5:350,1);
rho_V0 = z_rho_V0(200:5:350,2)/sum(z_rho_V0(200:5:350,2));% Important: always normalize rho here

% this index corresponds to the values of strain in the first column of the
% respective params files
%is = 1;

is_array = [1];%[1,2,3,4,9,14];

disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
for s = is_array
    strain = anisotropic_params(1, s)
end

%%%%%%%%%%%%%%%%
% function to generate all the parameters necessary for the
% calculations
% Takes a varying number of paramters, there are currently 2 numbers
% accepted:
% If you want to calculate graphene parameters for a minimum of
% V(0,0,z), the inputs are (strain_type, is(integer or array of intergers), "min")
% If you want to calculate graphene parameters for a specific value of
% z, the inputs are (strain_type, is(integer or array of intergers), "z", height(number>0))
% If you want to average the graphene potential in z 
% the inputs are (strain_type, is(integer or array of intergers), "average", [array of z's], [array of rhos])
%%%%%%%%%%%%%%%%

[V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
    Generate_graphene_parameters(strain_type, is_array, "average", z_V0, rho_V0);
%[V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = ...
%    Generate_graphene_parameters(strain_type, is_array, "z", 2);



