%  This function generates a number of the lowest shooting method solution
%  to the Schrodinger equation for some specified V. 
%  Inputs:
%  Psi0, dPsi0: guess at r=0
%  r: array of z's for which the solution is found
%  V: Potential
%  n: Number of solution to find
%  cutoff: how close to 0 we want the wavefunction at max r
%  initial_E, initial_step: initial guess for E and step
%  Outputs:
%  E_array: The energies found
%  Psi_found: The wavefunctions found
%  dPsi_found
function [E_array, Psi_found, dPsi_found] = shooting_method_function(Psi0, dPsi0, r, V, n, cutoff, initial_step, initial_E)
    foundn = 1; % Counts the number of solutions found
    E_find = initial_E; % Initial guess for energy
    step = initial_step; % Initial step for increasing energy
    zeroing = false; % Boolean for Bisection method

    % Containers for the found solutions
    Psi_found = zeros(n, length(V)); 
    dPsi_found = zeros(n, length(V)); 
    
    E_array = zeros(1,n); % The array of E solutions

    % Solve for the solution using the initial conditions
    [Psi, dPsi] = ShootingMethodPredict(Psi0, dPsi0, r, V, E_find);
    last_Psi = Psi(end);
    last_dPsi = dPsi(end);
    
    % Loop to find root
    while foundn<n+1 % While not converged
        E_find = E_find + step; % Increase E by the step
        [Psi, dPsi] = ShootingMethodPredict(Psi0, dPsi0, r, V, E_find); % Solve for Psi in the interval

        % Interval bisection method for finding E
        if (abs(Psi(end)) > cutoff)
            if not(zeroing) % Before the sign of the final step of Psi changes, just increase E
                if (last_Psi/Psi(end) < 0) % Check if the sign has changed
                    % Once the sign changes, immediately switches to the interval bisection step
                    zeroing = true;
                    step = step * -0.5;
                end
            else % Actual interval bisection steps
                if (last_Psi/Psi(end) < 0)
                    step = step * -0.5;
                else
                    step = step * 0.5;
                end
            end
                    
            last_Psi = Psi(end);
            last_dPsi = dPsi(end);
        else % Once the solution is found, get all outputs
            E_array(foundn) = E_find;
            
            Psi_found(foundn, :) = Psi;
            dPsi_found(foundn, :) = dPsi;
            
            E_find = E_find + initial_step;
            step = initial_step; % Reset the steps
            
            [Psi, dPsi] = ShootingMethodPredict(Psi0, dPsi0, r, V, E_find);
            last_Psi = Psi(end);
            last_dPsi = dPsi(end);
            
            foundn = foundn + 1;
            zeroing = false;
        end
    end
end