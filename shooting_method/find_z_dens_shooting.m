% Code to produce shooting method result for the He-Graphene potential

% load the He-Graphene potential in z
V_He_dat0 = dlmread("He_Graphene_V0.txt");
r_He_He0 = V_He_dat0(:,1); % First column is coordinates in z
V_He_He0 = V_He_dat0(:,2); % Second column is the potential

% Flip the potential
V_He_He_rev = flip(V_He_He0);
r_He_He_rev = flip(r_He_He0);

% Set up such that once the potential gets above max, it is just set to be
% the max
max = 0;
hit_max = false;
max_index = 0;

% The loop that sets the potential to the max once it hits that
for ir = 1:length(r_He_He0)
    if not(hit_max)
        if  V_He_He_rev(ir) > max
            hit_max = true;
            V_He_He_rev(ir) = max;
            max_index = ir;
        end
    else
    V_He_He_rev(ir) = max;
    end
end

% Truncates the potential given to the shooting potential code
z_cutoff = max_index + 50;

% Parameters for the shooting method function
psi0 = 0; % Initial Psi
dpsi0 = 1; % Initial dPsi
r = r_He_He0(1:z_cutoff); % range of z
V = V_He_He_rev(1:z_cutoff); % V(z)
n = 3; % Number of solutions to find
cutoff = 1; % Cutoff at the end of the interval for Psi
initial_step = 1; % Initial increase of E
initial_E = -250; % Initial guess for V

% Calls the shooting method function
[E_found, Psis, dPsis] = shooting_method_function(psi0, dpsi0, r, V, n, cutoff, initial_step, initial_E);

% Plots the wavefunctions
figure(1)
for in = 1:n
    subplot(n,1,in)
    plot(r_He_He_rev(1:z_cutoff), Psis(in,:))
    title("Unnormalized wavefunction")
    xlabel("z(Angstroms)")
end


% Put the wavefunction back on the original z array
for in = 1:n
    for ir = 1:z_cutoff
        Psis_extended(in, length(Psis_extended) + 1 - ir) = Psis(in,ir);
    end
end

rhos = 0*Psis_extended;

% Normalize the wavefunctions
for in = 1:n
    rhos(in,:) = (Psis_extended(in,:)).^2/sum((Psis_extended(in,:)).^2);
end

% Plots the densities
figure(2)
for in = 1:n
    subplot(n,1,in)
    plot(r_He_He0, rhos(in,:))
    title("Normalized Densities")
    xlabel("z(Angstroms)")
end

% Outputs the ground state Wavefunciton
dlmwrite("V0_rho0.txt", [r_He_He0 rhos(1,:)'])
