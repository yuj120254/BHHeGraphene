
% the paramters sigma and epsilon for all strain, using the latest values
params = [[0.00,0.05,0.10,0.15,0.16,0.17,0.18,0.19,0.20,0.21,0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.34],
          [2.643,2.669,2.700,2.738,2.749,2.759,2.769,2.780,2.791,2.804,2.817,2.831,2.846,2.858,2.878,2.896,2.917,2.938,3.065],
          [16.961,16.758,16.530,16.256,16.170,16.086,16.001,15.911,15.825,15.720,15.618,15.510,15.400,15.303,15.156,15.025,14.882,14.733,13.963]];
      
% loop for corrugation, currently unused
for i = 1:1%length(params[1])
    % sets the grid in z and C-C distance a0
    z = [2:0.01:6]; % In angstroms
    a0_real = 1.42; % In angstroms
    
    % strain
    delta = params(1, i);
    strain_string = sprintf("%02d", 100*params(1, i));
    
    poisson = 0.165; % poisson ratio, for anisotropic strain
    sigma = params(2, i); % In angstroms
    epsilon = params(3, i); % In Kelvin
    
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
    
    % generates V-He-Graphene at (x,y)=(0,0) for a range of z_nondim
    Va = V_z_0(z_nondim, 2, glat);

    % manually sets the z at which all the rest of the calculation happens
    zmin_nondim = 2.57*space_normalization;

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
    x = [-sqrt(3)*a0_nondim : sqrt(3)*a0_nondim/(mx) : 1.001*sqrt(3)*a0_nondim - sqrt(3)*a0_nondim/(mx)];
    y = [-1.5*a0_nondim : 1.5*a0_nondim/(my) : 1.001*1.5*a0_nondim - 1.5*a0_nondim/(my)];

    % calculates the He-Graphene potential in x and y
    Vr = real(V_fourier(x, y, gterms, Vfourier_scaled, 1, glat));

    V_dimensionless = Vr';

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
    MperiodsXlat = 12;
    MperiodsYlat = 12;

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

    % calculates potential in real units(angstroms) in cartesian space
    a0 = a0_real;
    x = 2*[-sqrt(3)*a0/2 : sqrt(3)*a0/(2*mx) : sqrt(3)*a0/2 - sqrt(3)*a0/(2*mx)];
    y = 2*[-1.5*a0/2 : 1.5*a0/(2*my) : 1.5*a0/2 - 1.5*a0/(2*my)];

    Vr = real(V_fourier(x, y, gterms, Vfourier, energy_scale*prefactor, glat));

    V_dimensionless = Vr';
    
    % writes output files
    dlmwrite(strcat("V_graphene", strain_string, ".txt"), V_dimensionless, ',')
    
    dlmwrite(strcat("V_graphene_lat", strain_string, ".txt"), V_dimensionless_lat, ',')
    
    dlmwrite(strcat("V_fourier", strain_string, ".txt"), real(Vfourier_scaled), ',')
        
    dlmwrite(strcat("Lattice_Vectors", strain_string, ".txt"), [glat.a1; glat.a2; glat.g1; glat.g2], ',')
    
    dlmwrite(strcat("Julia_grid", strain_string, ".txt"), [x; y]', ',')
    
    dlmwrite(strcat("Julia_params", strain_string, ".txt"), [space_normalization, energy_scale], ',')
end
