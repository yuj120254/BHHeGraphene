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

% this index corresponds to the values of strain in the first column of the
% respective params files
is = 1;

disp(strcat("computing graphene paramters for ", strain_type, " strain:"));
strain = anisotropic_params(1, is)

%%%%%%%%%%%%%%%%
% function to generate all the parameters necessary for the
% calculations
% Takes a varying number of paramters, there are currently 2 numbers
% accepted:
% If you want to calculate graphene parameters for a minimum of
% V(0,0,z), the inputs are (strain_type, is, "min")
% If you want to calculate graphene parameters for a specific value of
% z, the inputs are (strain_type, is, "z", 2.5), here 2.5 is the
% specified value, only change that.
%%%%%%%%%%%%%%%%
[V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = Generate_graphene_parameters(strain_type, is, "z", 2.57);

% writes output files
dlmwrite(strcat("V_graphene", strain_string, ".txt"), V, ',')

dlmwrite(strcat("V_graphene_lat", strain_string, ".txt"), V_lat, ',')

dlmwrite(strcat("V_fourier", strain_string, ".txt"), V_fc, ',')

dlmwrite(strcat("Lattice_Vectors", strain_string, ".txt"), lattice_vectors, ',')

dlmwrite(strcat("Julia_grid", strain_string, ".txt"), grid, ',')

dlmwrite(strcat("Julia_params", strain_string, ".txt"), param, ',')

% loop over different values of strain
for is = 1:1%length(params[1])

end
