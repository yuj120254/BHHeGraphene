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

function [V, V_lat, V_fc, lattice_vectors, grid, param, strain_string] = Graphene_paramters_function(varargin)
    % sets the grid in z 
    z = [2:0.01:6]; % In angstroms. This is only for minimizing V in z, not used in other calculations
    
    % poisson ratio
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
    iis = 1;% Index of strain array, useful when non consecutive strain values are used
    
    % loop over strain
    for is = is_array
        % establish parameters for isotropic of anisotropic strain
        if varargin{1} == "isotropic" 
            % C-C distance, affected by strain for isotropic strain
            a0_real = (1 + isotropic_params(1, is)) * 1.42; % In angstroms

            % only for anisotropic strain, does not change for isotropic strain
            delta = 0;
            strain_string = sprintf("i%02d", 100*isotropic_params(1, is)); % String for outputs

            sigma = isotropic_params(2, is); % In angstroms
            epsilon = isotropic_params(3, is); % In Kelvin
        elseif varargin{1} == "anisotropic" 
            % C-C distance
            a0_real = 1.42; % In angstroms

            % anisotropic strain
            delta = anisotropic_params(1, is);
            strain_string = sprintf("a%02d", 100*anisotropic_params(1, is)); % String for outputs

            sigma = anisotropic_params(2, is); % In angstroms
            epsilon = anisotropic_params(3, is); % In Kelvin
        end

        % lattice vectors
        a1_real = [(sqrt(3)*a0_real/2)*(1 - delta*poisson), (3*a0_real/2)*(1 + delta)]; % In angstroms
        a2_real = [-a1_real(1), a1_real(2)]; % In angstroms
        a3_real = a1_real-a2_real; % In angstroms

        % normalization for all lattice vectors such that their length is pi
        space_normalization = pi/norm(a1_real);
        z_nondim = z*space_normalization;
        sigma_nondim = sigma*space_normalization;
        a0_nondim = a0_real*space_normalization;
        a1_nondim = a1_real*space_normalization;
        a2_nondim = a2_real*space_normalization;
        a3_nondim = a3_real*space_normalization;
            
        % creates the Graphene_Lattice class for all subsequent calculations
        glat = Graphene_Lattice(sigma_nondim, delta, a0_nondim, a1_nondim, a2_nondim);

        % scaling factors for calculations
        prefactor = epsilon*glat.sigma*glat.sigma*2*pi/glat.A;

        h = 1.05457e-34; % kg m^2 /s (hbar)
        m = 6.6464767e-27; % kg

        E_R = h^2 * pi^2 / (2 * m * norm(a1_real*10^-10)^2); % k # ????24ma^2????
        k_B = 1.380649e-23; % J / K

        % scale for energy
        energy_scale = k_B/E_R;

        % chooses if we calculate the potential at the minimum of V in z,
        % or at the value of z input
        if nargin == 3
            % generates V-He-Graphene at (x,y)=(0,0) to find the minimum of z
            Va = V_z_0(z_nondim, 2, glat);

            [~, I] = min(Va);

            % z set to the minimum of V(0,0)
            zmin_nondim = z_nondim(I);

        elseif nargin == 4
            % manually sets the z at which all the rest of the calculation happens
            zmin_nondim = varargin{4}(iis)*space_normalization;
        else % for all other cases this is not important
            zmin_nondim = 2;
        end

        % sets the number of fourier components taken into account for the
        % calculations and calculates the fourier components
        gterms = 2;
        
        % Calculate the potential for anisotropic from an external
        % potential, taking the value of the SP and max of the potential
        if nargin == 6
            % Check to make sure the type of strain is correct
            if varargin{1} == "anisotropic" 
                disp("wrong type of strain")
                break
            end
            
            % set up the set of foureir coefficients that are the same
            Vf1_array = zeros(5,5);
            Vf1_array(2,2) = 1;
            Vf1_array(2,3) = 1;
            Vf1_array(3,2) = 1;
            Vf1_array(3,4) = 1;
            Vf1_array(4,3) = 1;
            Vf1_array(4,4) = 1;
            
            Vf2_array = zeros(5,5);
            Vf2_array(1,2) = 1;
            Vf2_array(2,1) = 1;
            Vf2_array(4,2) = 1;
            Vf2_array(2,4) = 1;
            Vf2_array(4,5) = 1;
            Vf2_array(5,4) = 1;

            % Set the locations at which the potential is sampled
            r_min = [0, 0];
            r_max = (a1_nondim+a2_nondim)/3;
            r_sp = a1_nondim/2;
            
            % Calculates the sum for if the Fourier coefficient is 1
            V1_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V1_max = real(V_fourier(r_max(1), r_max(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_max = real(V_fourier(r_max(1), r_max(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V1_sp = real(V_fourier(r_sp(1), r_sp(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_sp = real(V_fourier(r_sp(1), r_sp(2), gterms, Vf2_array, 1, glat))/energy_scale;
            
            % Set up the matrix  which goes from the coefficients to the
            % external potential
            a = V1_max - V1_min;
            b = V2_max - V2_min;
            c = V1_sp - V1_min;
            d = V2_sp - V2_min;

            A = [a b; c d]*energy_scale;
            
            % Sets up the values of the potential sampled
            external_potential_params = [varargin{6}(iis)-varargin{4}(iis), varargin{5}(iis)-varargin{4}(iis)]';
            
            % Calculates the Fourier coefficients
            Vf_fit = A\external_potential_params;
            
            % Constructs the full Fourier coefficients array
            Vfourier = (Vf_fit(1)*Vf1_array + Vf_fit(2)*Vf2_array);
            Vfourier_scaled = Vfourier*energy_scale;

        % Calculate the potential for isotropic from an external potential,
        % does not currently work
        elseif nargin == 8  
            % Check to make sure the type of strain is correct
            if varargin{1} == "isotropic" 
                disp("wrong type of strain")
                break
            end
            
            Vf1_array = zeros(5,5);
            Vf1_array(2,2) = 1;
            Vf1_array(4,4) = 1;
            
            Vf2_array = zeros(5,5);
            Vf2_array(2,3) = 1;
            Vf2_array(3,2) = 1;
            Vf2_array(3,4) = 1;
            Vf2_array(4,3) = 1;
            
            Vf3_array = zeros(5,5);
            Vf3_array(1,2) = 1;
            Vf3_array(2,1) = 1;
            Vf3_array(4,5) = 1;
            Vf3_array(5,4) = 1;

            Vf4_array = zeros(5,5);
            Vf4_array(4,2) = 1;
            Vf4_array(2,4) = 1;
            
            r_min = [0, 0];
            r_half_max_para = glat.b2;%(a1_nondim+a2_nondim)/3;
            r_half_sp_para = a1_nondim/2;
            r_half_max_perp = glat.b1;%a3_nondim/2 + r_max_para/2;
            r_half_sp_perp = a3_nondim/2;
%             r_half_max_para = a1_nondim/2;%(a1_nondim+a2_nondim)/3;
%             r_half_sp_para = a1_nondim/2;
%             r_half_max_perp = a3_nondim/4;%a3_nondim/2 + r_max_para/2;
%             r_half_sp_perp = a3_nondim/2;

            V1_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V3_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf3_array, 1, glat))/energy_scale;
            V4_min = real(V_fourier(r_min(1), r_min(2), gterms, Vf4_array, 1, glat))/energy_scale;
            V1_half_max_para = real(V_fourier(r_half_max_para(1), r_half_max_para(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_half_max_para = real(V_fourier(r_half_max_para(1), r_half_max_para(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V3_half_max_para = real(V_fourier(r_half_max_para(1), r_half_max_para(2), gterms, Vf3_array, 1, glat))/energy_scale;
            V4_half_max_para = real(V_fourier(r_half_max_para(1), r_half_max_para(2), gterms, Vf4_array, 1, glat))/energy_scale;
            V1_half_max_perp = real(V_fourier(r_half_max_perp(1), r_half_max_perp(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_half_max_perp = real(V_fourier(r_half_max_perp(1), r_half_max_perp(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V3_half_max_perp = real(V_fourier(r_half_max_perp(1), r_half_max_perp(2), gterms, Vf3_array, 1, glat))/energy_scale;
            V4_half_max_perp = real(V_fourier(r_half_max_perp(1), r_half_max_perp(2), gterms, Vf4_array, 1, glat))/energy_scale;
            V1_half_sp_para = real(V_fourier(r_half_sp_para(1), r_half_sp_para(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_half_sp_para = real(V_fourier(r_half_sp_para(1), r_half_sp_para(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V3_half_sp_para = real(V_fourier(r_half_sp_para(1), r_half_sp_para(2), gterms, Vf3_array, 1, glat))/energy_scale;
            V4_half_sp_para = real(V_fourier(r_half_sp_para(1), r_half_sp_para(2), gterms, Vf4_array, 1, glat))/energy_scale;
            V1_half_sp_perp = real(V_fourier(r_half_sp_perp(1), r_half_sp_perp(2), gterms, Vf1_array, 1, glat))/energy_scale;
            V2_half_sp_perp = real(V_fourier(r_half_sp_perp(1), r_half_sp_perp(2), gterms, Vf2_array, 1, glat))/energy_scale;
            V3_half_sp_perp = real(V_fourier(r_half_sp_perp(1), r_half_sp_perp(2), gterms, Vf3_array, 1, glat))/energy_scale;
            V4_half_sp_perp = real(V_fourier(r_half_sp_perp(1), r_half_sp_perp(2), gterms, Vf4_array, 1, glat))/energy_scale;

            a = V1_half_sp_perp - V1_min;
            b = V2_half_sp_perp - V2_min;
            c = V3_half_sp_perp - V3_min;
            d = V4_half_sp_perp - V4_min;
            e = V1_half_sp_para - V1_min;
            f = V2_half_sp_para - V2_min;
            g = V3_half_sp_para - V3_min;
            h = V4_half_sp_para - V4_min;
            i = V1_half_max_perp - V1_min;
            j = V2_half_max_perp - V2_min;  
            k = V3_half_max_perp - V3_min;
            l = V4_half_max_perp - V4_min;
            m = V1_half_max_para - V1_min;
            n = V2_half_max_para - V2_min;
            o = V3_half_max_para - V3_min;
            p = V4_half_max_para - V4_min;

            A = [a b c d; e f g h; i j k l; m n o p]*energy_scale
            
            external_potential_params = [varargin{5}(iis)-varargin{4}(iis), varargin{6}(iis)-varargin{4}(iis),...
                                         varargin{7}(iis)-varargin{4}(iis), varargin{8}(iis)-varargin{4}(iis)]'
                                                 
            Vf_fit = A\external_potential_params
            
            Vfourier = (Vf_fit(1)*Vf1_array + Vf_fit(2)*Vf2_array + ...
                        Vf_fit(3)*Vf3_array + Vf_fit(4)*Vf4_array);
            Vfourier_scaled = Vfourier*energy_scale
        else
            Vfourier = Fourier_Terms(gterms, zmin_nondim, glat);
            Vfourier_scaled = Vfourier*energy_scale*prefactor;
        end

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

        xlat_x_array = [-Lxlat_x/2 : dxlat_x : Lxlat_x/2 + dxlat_x/4];    
        xlat_y_array = [-Lxlat_y/2 : dxlat_y : Lxlat_y/2 + dxlat_y/4]; 
        ylat_x_array = [-Lylat_x/2 : dylat_x : Lylat_x/2 + dylat_x/4]; 
        ylat_y_array = [-Lylat_y/2 : dylat_y : Lylat_y/2 + dylat_y/4]; 

        Mx_ones = ones(length(xlat_x_array));
        My_ones = ones(length(ylat_x_array));

        xgrid = zeros(length(xlat_x_array),length(xlat_x_array));
        ygrid = zeros(length(ylat_x_array),length(ylat_x_array));

        for ix = 1:size(xgrid, 1)
            for iy = 1:size(ygrid, 2)
                xgrid(ix, iy) = xlat_x_array(ix) + ylat_x_array(iy);
                ygrid(ix, iy) = xlat_y_array(ix) + ylat_y_array(iy);
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
            zs = varargin{4}(iis,:)*space_normalization; % make sure to convert to dimensionless coords
            rhos = varargin{5}(iis,:);

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

        figure(is)
        surf(x, y, Vr')
        xlabel('x')
        ylabel('y')
        view(0,90)

        % writes output files
        dlmwrite(strcat("V_graphene", strain_string, ".txt"), V, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("V_graphene_lat", strain_string, ".txt"), V_lat, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("V_fourier", strain_string, ".txt"), real(V_fc), 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Lattice_Vectors", strain_string, ".txt"), lattice_vectors, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Julia_grid", strain_string, ".txt"), grid, 'delimiter', ',', 'precision', 16)

        dlmwrite(strcat("Julia_params", strain_string, ".txt"), param, 'delimiter', ',', 'precision', 16)
        
        iis = iis + 1;
    end
end