# BHHeGraphene

Data generation for the single particle sections of: https://arxiv.org/abs/2102.11288

Currently, three files are needed to run the main file for the graphene parameters generation:
Graphene_lattice.m contains the class necessary for the calculations
Graphene_parameters_function.m contains the function used to generate the graphene parameters and output them to file.
V0_rho0.txt contains an example of a He density in z.
None of these files should be modified.

All of the important function calls and parameter setting is in the file: Generate_graphene_parameters_main.m. This should be the only file you modify to generate different graphene parameters. All of the instructions necessary to run the file is contained within the file itself.
Currently there are 4 working cases of He on Graphene that is included in this file:
1) Calculating Graphene Parameters at the minimum of V(0,0,z) in z
2) Calculating Graphene Parameters at some fixed z
3) Calculating Graphene Parameters while averaging in z
4) Calculating Graphene Parameters for some external potential with isotropic strain.


Currently the file make_2D_Bloch_bands_v13a.m is set up to work with the example of case 1 in the main file.

The files in the shooting method folder produces the ground state wavefunction for the third case of the main file.
The only file that needs to be ran in that file is find_z_dens_shooting.m, but all of the .m files are required to run it.
It takes some input in the form of two columns, z in Angstroms and V(z). The example given is the file called He_Graphene_V0.txt.
It outputs the ground state wavefunction of the solution to that potential, in two columns: z in Angstroms and rho(z), whose name can be changed. 
The example given is V0_rho0.txt. This file is then to be copied into the same folder as Graphene_parameters_function.m when running the 3rd case.
For each value of strain, the density file must be individually loaded and normalized, then the densities are put into the array as explained in the file. 
