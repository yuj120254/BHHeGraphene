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
n = max_index + 50;
figure(1)
plot(r_He_He0(1:n), V_He_He_rev(1:n));

% Calls the shooting method function
[E_found, Psis, dPsis] = shooting_method_function(0, 1, r_He_He0(1:n), V_He_He_rev(1:n), 1, 1, 1, -250);

% Plots the first wavefunction
figure(3)
plot(r_He_He_rev(1:n), Psis(1,:))

% Put the wavefunction back on the original z array
Psi1 = 0*r_He_He0;
for ir = 1:n
    Psi1(length(Psi1) + 1 - ir) = Psis(1,ir);
end

% Outputs the Wavefunciton
dlmwrite("V0_rho0.txt", [r_He_He0 Psi1])
