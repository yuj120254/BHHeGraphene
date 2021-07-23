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

function [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = Generate_graphene_paramters(varargin)
    % sets the grid in z and C-C distance a0
    z = [2:0.01:6]; % In angstroms This is only for minimizing V in z, not used in other options
    a0_real = 1.42; % In angstroms
    poisson = 0.165; % poisson ratio, for anisotropic strain
    
    % the paramters sigma and epsilon for all strain, using the latest values
    anisotropic_params = [[0.00,0.05,0.10,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.34],
                      [2.643,2.669,2.700,2.738,2.749,2.759,2.769,2.780,2.791,2.804,2.817,2.831,2.846,2.858,2.878,2.896,2.917,2.938,3.065],
                      [16.961,16.758,16.530,16.256,16.170,16.086,16.001,15.911,15.825,15.720,15.618,15.510,15.400,15.303,15.156,15.025,14.882,14.733,13.963]];
 
    isotropic_params = [[0.00,0.05,0.10,0.15,0.20,0.25,0.30],
                      [2.642,2.697,2.749,2.800,2.849,2.897,2.942],
                      [16.968,16.585,16.248,15.952,15.693,15.463,15.259]];
   
    % array of strain values that we are calculating parameters for
    is_array = varargin{2};
    
    % loop over strain
    for is = is_array
        % establish parameters for isotropic of anisotropic strain
        if varargin{1} == "isotropic" 
            % C-C distance
            a0_real = (1 + isotropic_params(1, is)) * 1.42; % In angstroms

            % strain
            delta = 0;
            strain_string = sprintf("%02d", 100*isotropic_params(1, is));

            sigma = isotropic_params(2, is); % In angstroms
            epsilon = isotropic_params(3, is); % In Kelvin
        elseif varargin{1} == "anisotropic" 
            % C-C distance
            a0_real = 1.42; % In angstroms

            % strain
            delta = anisotropic_params(1, is);
            strain_string = sprintf("%02d", 100*anisotropic_params(1, is));

            sigma = anisotropic_params(2, is); % In angstroms
            epsilon = anisotropic_params(3, is); % In Kelvin
        end

        % lattice vectors
        a1_real = [(sqrt(3)*a0_real/2)*(1 - delta*poisson), (3*a0_real/2)*(1 + delta)]; % In angstroms
        a2_real = [-a1_real(1), a1_real(2)]; % In angstroms

        % normalization for all lattice vectors such that their length is pi
        space_normalization = pi/norm(a1_real);
        z_nondim = z*space_normalization;
        sigma_nondim = sigma*space_normalization;
        a0_nondim = a0_real*space_normalization;
        a1_nondim = a1_real*space_normalization;
        a2_nondim = a2_real*space_normalization;

        % creates the Graphene_Lattice class for all subsequent calculations
        glat = Graphene_Lattice(sigma_nondim, delta, a0_nondim, a1_nondim, a2_nondim);

        % scaling factors for calculations
        prefactor = epsilon*glat.sigma*glat.sigma*2*pi/glat.A;

        h = 1.05457e-34; % kg m^2 /s (hbar)
        m = 6.6464767e-27; % kg

        E_R = h^2 * pi^2 / (2 * m * norm(a1_real*10^-10)^2); % k # ????24ma^2????
        k_B = 1.380649e-23; % J / K

        % chooses if we calculate the potential at the minimum of V in z,
        % or at the value of z input
        if nargin == 3
            % generates V-He-Graphene at (x,y)=(0,0) to find the minimum of z
            Va = V_z_0(z_nondim, 2, glat);

            [~, I] = min(Va);

            zmin_nondim = z_nondim(I);

        elseif nargin == 4
            % manually sets the z at which all the rest of the calculation happens
            zmin_nondim = varargin{4}*space_normalization;
        else % just some random value so the next part of the code runs
            zmin_nondim = 2;
        end

        % sets the number of fourier components taken into account for the
        % calculations and calculates the fourier components
        gterms = 2;
        Vfourier = Fourier_Terms(gterms, zmin_nondim, glat);

        % scale the fourier components
        energy_scale = k_B/E_R;
        Vfourier_scaled = Vfourier*energy_scale*prefactor;

        % sets how fine the grid is for each unit cell
        mx = 24;
        my = 24;
        % mulitplied by 1.001 to make sure the range reaches the end
        % x = [-sqrt(3)*a0_nondim : sqrt(3)*a0_nondim/(mx) : 1.001*sqrt(3)*a0_nondim - sqrt(3)*a0_nondim/(mx)];
        % y = [-1.5*a0_nondim : 1.5*a0_nondim/(my) : 1.001*1.5*a0_nondim - 1.5*a0_nondim/(my)];

        % calculates the He-Graphene potential in x and y
        % Vr = real(V_fourier(x, y, gterms, Vfourier_scaled, 1, glat));

        % V_dimensionless = Vr';

        % change of variables to basis with axis aligned along the lattice
        % vectors
        xlat_x = glat.a1(1); 
        xlat_y = glat.a1(2); 
        ylat_x = glat.a2(1);
        ylat_y = glat.a2(2);

        % sets how fine the grid is along each axis
        Mx_singlePeriod = 12;
        My_singlePeriod = 12;

        % calculates the step along each axis
        dxlat_x = xlat_x / Mx_singlePeriod; 
        dxlat_y = xlat_y / Mx_singlePeriod; 
        dylat_x = ylat_x / My_singlePeriod;
        dylat_y = ylat_y / My_singlePeriod;

        % number of unit cells along each axis
        % equal to 2x Qx Qy in make_bloch_bands
        MperiodsXlat = 24;
        MperiodsYlat = 24;

        % makes grid in new basis
        Lxlat_x = MperiodsXlat * xlat_x; 
        Lxlat_y = MperiodsXlat * xlat_y; 
        Lylat_x = MperiodsYlat * ylat_x;
        Lylat_y = MperiodsYlat * ylat_y;

        xlat_x = [-Lxlat_x/2 : dxlat_x : Lxlat_x/2 + dxlat_x/4];    
        xlat_y = [-Lxlat_y/2 : dxlat_y : Lxlat_y/2 + dxlat_y/4]; 
        ylat_x = [-Lylat_x/2 : dylat_x : Lylat_x/2 + dylat_x/4]; 
        ylat_y = [-Lylat_y/2 : dylat_y : Lylat_y/2 + dylat_y/4]; 

        Mx_ones = ones(length(xlat_x));
        My_ones = ones(length(ylat_x));

        xgrid = zeros(length(xlat_x),length(xlat_x));
        ygrid = zeros(length(ylat_x),length(ylat_x));

        for ix = 1:size(xgrid, 1)
            for iy = 1:size(ygrid, 2)
                xgrid(ix, iy) = xlat_x(ix) + ylat_x(iy);
                ygrid(ix, iy) = xlat_y(ix) + ylat_y(iy);
            end
        end

        % calculates the potential in the xlat-ylat grid
        Vr_lat = real(V_fourier_grid(xgrid, ygrid, gterms, Vfourier_scaled, 1, glat));

        V_dimensionless_lat = Vr_lat';

        a0 = a0_real;
        x = 2*[-sqrt(3)*a0/2 : sqrt(3)*a0/(2*mx) : sqrt(3)*a0/2 - sqrt(3)*a0/(4*mx)];
        y = 2*[-1.5*a0/2 : 1.5*a0/(2*my) : 1.5*a0/2 - 1.5*a0/(4*my)];

        Vr = real(V_fourier(x, y, gterms, Vfourier, energy_scale*prefactor, glat));

        V_dimensionless = Vr';

        % for the averaging in z case
        if nargin == 5
            % sets the arrays where the potential and fc's are accumulated 
            V = 0*V_dimensionless;
            V_lat = 0*V_dimensionless_lat;
            V_fc = 0*real(Vfourier_scaled);

            % gets the z and rho
            zs = varargin{4}*space_normalization; % make sure to convert to dimensionless coords
            rhos = varargin{5};

            % loops over z
            for iz = 1:length(zs)
                
                % sets the number of fourier components taken into account for the
                % calculations and calculates the fourier components
                gterms = 2;
                Vfourier = Fourier_Terms(gterms, zs(iz), glat);

                % scale the fourier components
                energy_scale = k_B/E_R;
                Vfourier_scaled = Vfourier*energy_scale*prefactor;
                V_fc = V_fc + real(Vfourier_scaled)*rhos(iz); % adds to weighted average

                % calculates the potential in the xlat-ylat grid
                Vr_lat = real(V_fourier_grid(xgrid, ygrid, gterms, Vfourier_scaled, 1, glat));

                V_dimensionless_lat = Vr_lat';
                V_lat = V_lat + V_dimensionless_lat*rhos(iz); % adds to weighted average

                % calculate the potential on a square grid
                Vr = real(V_fourier(x, y, gterms, Vfourier, energy_scale*prefactor, glat));

                V_dimensionless = Vr';
                V =  V + V_dimensionless*rhos(iz); % adds to weighted average
            end
        else
            % the outputs when not averaging in z
            V = V_dimensionless;
            V_lat = V_dimensionless_lat;
            V_fc = real(Vfourier_scaled);
        end
        
        % outputs
        lattice_vectors = [glat.a1; glat.a2; glat.g1; glat.g2];
        grid = [x; y];
        param = [space_normalization, energy_scale];

        % writes output files
        dlmwrite(strcat("V_graphene", strain_string, ".txt"), V, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("V_graphene_lat", strain_string, ".txt"), V_lat, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("V_fourier", strain_string, ".txt"), real(V_fc), 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Lattice_Vectors", strain_string, ".txt"), lattice_vectors, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Julia_grid", strain_string, ".txt"), grid, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Julia_params", strain_string, ".txt"), param, 'delimiter', ',', 'precision', 16)
    end
end