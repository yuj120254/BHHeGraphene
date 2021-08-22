%  This function generates a shooting method solution to the Schrodinger 
%  equation for some specified V for a specific E.
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
function [x1, x2] = ShootingMethodPredict(Psi0, dPsi0, r, V, E)
    % Initialize the solution and its derivative
    x1 = zeros(1, length(V));
    x2 = zeros(1, length(V));
    x1(1) = Psi0;
    x2(1) = dPsi0;
    
    % The step
    h = r(2) - r(1);
        
    % Interpolates a V on a 2x finer grid
    r_finer = r(1):h/2:r(end);
    V_finer = spline(r, V, r_finer);
    
    % Runge-Kutta method
    for i = 1:length(V)-1
        k11 = x2(i);
        k21 = f2(x1(i), V_finer(2*i), E);
        k12 = x2(i) + h*k21/2;
        k22 = f2(x1(i) + h*k11/2, V_finer(2*i), E);
        k13 = x2(i) + h*k22/2;
        k23 = f2(x1(i) + h*k12/2, V_finer(2*i), E);
        k14 = x2(i) + h*k23;
        k24 = f2(x1(i) + h*k13, V_finer(2*i), E);
        x1(i+1) = x1(i) + h*(k11 + 2*k12 + 2*k13 + k14)/6;
        x2(i+1) = x2(i) + h*(k21 + 2*k22 + 2*k23 + k24)/6;
    end
end