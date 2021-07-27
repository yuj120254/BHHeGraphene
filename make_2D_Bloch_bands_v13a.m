% This code determines 2D bant)dgap structure in Fig. 7.5 for Eq. (7.39): 
% Uxx + Uyy + mu*U + V(x,y)U = 0.
% It then computes the Wannier function corresponding to a given energy band.

% Changelog for version 4, 5, 6, and 6.5 is at the very end of this code. 

clear all

strains = [0.00,0.05,0.10,0.15,0.16,0.17,0.18,0.19,0.20,0.21,...
           0.22,0.23,0.24,0.25,0.26,0.27,0.28,0.29,0.34];
       
strain_strings = ["a00", "a05", "a10", "a15", "a16", "a17", "a18", "a19", "a20", "a21",...
                  "a22", "a23", "a24", "a25", "a26", "a27", "a28", "a29", "a34"]; 

% strains = [0.00,0.05,0.10,0.15,0.20,0.25,0.29,0.34];
%        
% strain_strings = ["00", "05", "10", "15", "20", "25", "29", "34"];

strain_index_array = [1,2,3,4,9,14];
           
for iss = 1:1%strain_index_array
    %for iss = 1:1%length(strains)
    compare_w1D = 0;  % 1  =>  Compare our 2D results with Adrian's 1D results.
                      % 0  =>  Do not compare. 

    % Load Params
    Julia_Params = dlmread(strcat('Julia_params', ...
                                  strain_strings(iss),'.txt'));
    space_normalization = Julia_Params(1);
    energy_normalization = Julia_Params(2);

    cutoff = 1.92*space_normalization;

    % Load the g vectors
    Lattice_Vectors = dlmread(strcat('Lattice_Vectors', ...
                                     strain_strings(iss) ,'.txt'));

    %a0 = 1.42;
    a1 = Lattice_Vectors(1,:);
    a2 = Lattice_Vectors(2,:);
    g1 = Lattice_Vectors(3,:);
    g2 = Lattice_Vectors(4,:);

    a1_nondim = a1;
    a2_nondim = a2;
    g1_nondim = g1;
    g2_nondim = g2;

    % Nondimensional areas of the Brillouin Zone and the Primary Unit Cell:
    A_BZ = abs(det( [g1_nondim.'  g2_nondim.'] ));
    A_PUC = abs(det( [a1_nondim.'  a2_nondim.'] ));
    AA_check = A_BZ * A_PUC / (2*pi)^2;
    disp(['ABZ = ' num2str(A_BZ) ';  APUC = ' num2str(A_PUC) '; ABZ*APUC/(2pi)^2 = ' num2str(AA_check)])
    if  abs(AA_check - 1) > 1e-5
        disp('Stop: ABZ * APUC ~= (2pi)^2;  There is an ERROR somewhere !')
        adsf
    end

    a1norm = a1/norm(a1);
    a2norm = a2/norm(a2);
    a1a2norm_dot = a1norm*a2norm';

    norm_g1_nondim = norm(g1_nondim);
    norm_g2_nondim = norm(g2_nondim);
    g1norm = g1/norm_g1_nondim;
    g2norm = g2/norm_g2_nondim;
    cosine_g1g2 = (g1_nondim * g2_nondim')/(norm_g1_nondim * norm_g2_nondim);

    % --- 
    % Set up the potential in a single period:
    Lxlat_singlePeriod = pi; 
    Lylat_singlePeriod = pi*norm(a2)/norm(a1);  % x- and y-lengths of the cell
    Mxlat_singlePeriod = 12;  % xlat || a1, ylat || a2
    Mylat_singlePeriod = 12;  
    % !!!!!Updating these requires updating them in julia!!!!!!
                           % # of discretization points per 1 period 
                           % For technical reasons, take them to be *even*.
                           % Then 0 occurs at a grid point. 

    dxlat_singlePeriod = Lxlat_singlePeriod / Mxlat_singlePeriod; 
    dylat_singlePeriod = Lylat_singlePeriod / Mylat_singlePeriod;
    xlat_singlePeriod = -Lxlat_singlePeriod/2 : dxlat_singlePeriod : Lxlat_singlePeriod/2 - dxlat_singlePeriod; 
    ylat_singlePeriod = -Lylat_singlePeriod/2 : dylat_singlePeriod : Lylat_singlePeriod/2 - dylat_singlePeriod;
                   % For purely technical reasons, related to the subsequent calculation of 
                   %  Wannier functions, we define the (x,y)-cell to be (almost) centered at 0. 
    [Xlat_singlePeriod,  Ylat_singlePeriod] = meshgrid(xlat_singlePeriod, ylat_singlePeriod); 

    % Loads the pre calculated fourier coefficents of V
    Vfourier = dlmread(strcat('V_fourier', strain_strings(iss) ,'.txt'));
    V_multiplier = 1;
    V_sign = -1;
    Vfourier = V_sign*V_multiplier*Vfourier;
    Vf_size = size(Vfourier);
    gterms = (Vf_size(1)-1)/2;

    % "pads" the Vfourier with the extra zeros to get the required number of
    % fourier coefficents
    increase_factor = 1;
    % For the strength of the C-He potential that is appropriate for graphene,
    % using increase-factor = 1 is adequate. This can be seen from the log plot
    % of the Wannier function, which has the correct 6-fold symmetry.
    % However, as the potential's strength (= depth) increases to be more than
    % 50% the correct one, the Wannier function develops a 4-fold symmetry
    % far from the center. This indicates that we need more Fourier coefficients
    % to compute it correctly.
    % For the strength = 2 x the "regular" potential, using  increase_factor = 2
    % instead of  increase_factor = 1  restores the correct 6-fold symmetry of the WF.
    Vfourier_big = zeros(increase_factor*4*gterms + 1, increase_factor*4*gterms + 1);
    for i = 1 : 2*gterms + 1
        for j = 1 : 2*gterms + 1
            Vfourier_big(i + (2*increase_factor-1)*gterms, ...
                         j + (2*increase_factor-1)*gterms) = Vfourier(i,j);
        end
    end
    Vfourier = Vfourier_big;

    % For debugging only
    % Vfourier;
    % 
    % pause



    % --- 
    % Set up the Brillouin zone:
    %           VERY IMPORTANT NOTE:
    % - For the computation of energy vs wavenumber, we must use limits |g{1,2}_nondim| (divided by 2),
    %   because this is what appears in the equation for the Fourier coefficients and
    %   this is also the side of the parallelogram representing the Brilloin zone.
    % - However, in the situations where we have to compute terms like 
    %   exp( i * vec{k} \dot \vec{r} ),
    %   we must set DIFFERENT limits for \vec{g} (i.e., "k0x" and "k0y") in order to ENSURE
    %   that max| \vec{k} \dot \vec{r} | = \pi. 
    %   The issue stems from the fact that \vec{g_i} is NOT orthogonal to \vec{a_i}
    %   and hence |g_i| * |a_i| > 2*pi. 
    %   Thus, when we compute terms like  exp( i * vec{k} \dot \vec{r} ),
    %   we must use new_k0{x,y} = 2*pi/L{x,y}lat_singlePeriod .
    % - So:
    %      1) We will define k-vectors to be on the intervals [ -k0{x,y}/2, k0{x,y}/2 ),
    %   where  k0{x,y} = 2*pi/L{x,y}lat_singlePeriod,  to guarantee condition  
    %   max| \vec{k} \dot \vec{r} | = \pi.
    %   These will be used ONLY in the aforementioned exp-terms.
    %      2) We will define a "non-rectangularity factor":
    %   nonrectang_{x,y} = |g{1,2}_nondim| / k0{x,y}.
    %   The k-values IN THE ENERGY CALCULATIONS ONLY will need to be multiplied by it in order
    %   to have k-values run in the range [ -|g{1,2}_nondim|/2, |g{1,2}_nondim| ) 
    %   instead of in the range defined in 1). 
    %
    %   1) "Length" and grids in the reciprocal (=wavenumber) space PROJECTED on the 
    % direct lattice vectors (i.e. on {a1, a2}) and used only in the exp-like terms (see above).
    % % % % Lengths and grids along the sides of the BZ parallelogram,
    % % % % used in the energy(k) calculations (see above):
    k0x = 2*pi / Lxlat_singlePeriod; 
    k0y = 2*pi / Lylat_singlePeriod; % width of Brillouin zones

    Qx = 12; 
    Qy = 12;   % 1/2*(number of grid points along edges of Brillouin zone(MperiodsX/Ylat in Julia))  
    % kxgrid = -(k0x/2)+ 0.5*(k0x/2)/Qx: (k0x/2)/Qx : (k0x/2) - 0.5*(k0x/2)/Qx; 
    % kygrid = -(k0y/2) + 0.5*(k0y/2)/Qy : (k0y/2)/Qy : (k0y/2) - 0.5*(k0y/2)/Qy;
    %%%%%%%%%%%%%%
    %%% - 10c ---> 
    dkxgrid = k0x/(2*Qx);
    dkygrid = k0y/(2*Qy);
    kxgrid = ( -(k0x/2) : dkxgrid : (k0x/2) - dkxgrid )     + dkxgrid; 
    kygrid = ( -(k0y/2) : dkygrid : (k0y/2) - dkygrid )     + dkygrid;
    %
    %     2) Define the "nonrectangularity factors" defined in 2) above;
    %  they are needed to obtain correct lengths (i.e., |g{1,2}_nondim|) along the sides 
    %  of the BZ parallelogram, and used ONLY in the energy(k) calculations.
    nonrectang_x = norm_g1_nondim / k0x; 
    nonrectang_y = norm_g2_nondim / k0y;


    %%%% <--- 10c -
    %%%%%%%%%%%%%%
    % % The issue addressed in the Note below has not been seen in any cases where V \neq 0.
    % % In fact, we now DO take k-values on the edge of the Brilloin zone. 
    %        % Explanation:
    %        % We avoid taking points on the edge of the Brillouin zone, where kx or ky = pi/L. 
    %        % Reason:
    %        % Later we compute eigenvalues of a matrix, called  HamiltonianMatrix  below, 
    %        % that involves terms  (kx + Km)^2 and (ky + Kn)^2, 
    %        % where {Km,Kn} = {m,n}*2*pi/{Lx,Ly}, {m,n}=integer; and "single period" for {Lx,Ly}
    %        % is implied everywhere in this note. 
    %        % Then for kx = pi/L_x,  *two values* of Km will give the same (kx + Km)^2:
    %        %  Km1 = 0 and Km2 = -2*pi/L_x. 
    %        % Then if the matrix elements of V also happen to have some symmetry involving
    %        % indices m1 = 0 and m2 = -1 (as, e.g., happens in the trivial case V = 0),
    %        % then  HamiltonianMatrix  will have repeated eigenvalues, which will then have
    %        % different eigenvectors for the same eigenvalue. This will create a technical
    %        % difficulty in selecting eigenvectors that *continuously* depend on {kx,ky}
    %        % (as we want only 1 eigenvector per energy band). 
    %        % Thus, by selecting kx, ky that never are integer*pi/{Lx,Ly}, we avoid having 
    %        % double eigenvalues. 
           %


    % ---
    % Set up auxiliary parameters related to the computation of the eigenvalues
    %  and the Bloch functions:
    Nx = increase_factor * 2*gterms;  
         % increase-factor  was defined before  Vfourier_big  some 50-100 lines above
    Ny = increase_factor * 2*gterms;    % (2Nx+1) and (2Ny+1) are the numbers of Fourier modes along x and y 
               %  in expansion of the periodic potential V and the corresponding 
               %  periodic part, G, of the Bloch functions.


    Kmvec = k0x*(-Nx : Nx);
    Knvec = k0y*(-Ny : Ny);   % Fourier k-vectors needed for these calculations. 
    % % % [Km, Kn] = meshgrid(Kmvec, Knvec);
    [Kn, Km] = meshgrid(Knvec, Kmvec);   % matrix of dimension (2Nx+1) x (2Ny+1) 


    % --- 
    % Set up (x,y)-domain that includes several periods, to compute Wannier functions:
    MperiodsXlat = 2*Qx; 
    MperiodsYlat = 2*Qy;
                 % For the computatin of Wannier functions (which are defined on the
                 % infinite interval), the number of periods per dimension must equal
                 % the number of k-values per one side of the Brilloin zone. 
                 % (For the energy calculations alone, no space-related calculations 
                 %  are needed, and so the number of periods simply doesn't appear.)
    Lxlat = MperiodsXlat * Lxlat_singlePeriod; 
    Lylat = MperiodsYlat * Lylat_singlePeriod;   % length of x- and y-sides
    xlat = -Lxlat/2 : dxlat_singlePeriod : Lxlat/2 - 1*dxlat_singlePeriod;    
    ylat = -Lylat/2 : dylat_singlePeriod : Lylat/2 - 1*dylat_singlePeriod; 
       %  The setup with "0*..." instead of "1*..." in front of the last term,
       %  which ensures that  min(x,y) = -max(x,y),  whereby one has the same number 
       %  of grid points on both sides of 0. 
       %  This will be convenient for Wannier function (WF) calculations if in a later version
       %  of the code we decide to use the symmetry of the potential to compute WF using 
       %  the pseudo-momenta {kx,ky} from only the 1st quadrant of the Brillouin zone. 
    Mxlat = length(xlat);
    Mylat = length(ylat);   % number of grid points in the multi-period x- and y-domains


    % ---
    % Center of Wanner function:
    siteX = 0; 
    siteY = 0;  


    % ---
    % Set up exponentials that will be used to calculate Wannier functions out of Bloch ones:
    for ix = 1:Mxlat
        eiKmx(:,:, ix) = exp(1i*Km*xlat(ix));
    end
    %
    for iy = 1:Mylat
        eiKny(:,:, iy) = exp(1i*Kn*ylat(iy));
    end
    % Both of these are matrixes of dimensions (2Nx+1) x (2Ny+1). 


    % ---
    % Calculate Fourier coefficients of potential V (see the NOTE after this loop):
    % Vfourier = zeros(2*Nx+1, 2*Ny+1);
    % for m = 1 : 2*Nx+1 
    %   for n = 1 : 2*Ny+1 
    % %     [m n]
    %     integrand = V.*exp( -i*( (m-Nx-1)*k0x*X_singlePeriod + (n-Ny-1)*k0y*Y_singlePeriod ) ); 
    %     Vfourier(m, n) = dx_singlePeriod*dy_singlePeriod * sum(sum(integrand))/...
    %                     (Lx_singlePeriod*Ly_singlePeriod);
    % % % %     Vfourier(n, m) = dx_singlePeriod*dy_singlePeriod * sum(sum(integrand))/...
    % % % %                     (Lx_singlePeriod*Ly_singlePeriod);
    % 
    % %     ddd = Vfourier(m,n)   % for debugging only
    % %     pause
    %   end
    % end 
    %         % This is for debugging only:
    %         aaa = real(Vfourier)
    %
    %%%
    % NOTE: 
    % In a later edition of this code we may consider using twice as many (per dimension)
    % Fourier coefficients of V than we will use coefficients of G (the periodic part
    % of the Bloch function). 
    % For now this doesn't matter too much since Fourier coefficients of smooth potentials,
    % that we have used so far, decay very rapidly with their Fourier index (i.e., wavenumber).

    % For debugging only: 
    % figure(1000)
    % mesh(Kmvec/2,Knvec/2, abs(Vfourier));
    % xlabel('Km'); ylabel('Kn'); zlabel('|Vfourier|')



    disp('Fourier coefficients of V have been computed.')
    disp('Working on finding the matrix form of the original equation.')

    tic

    % -------------------------------
    % Begin calculation of the energy bands, after which compute Wannier functions.

    Num_mu = 4;  % number of energy bands that we will seek

    % Computing the matrix form of the Hamiltonian and from it,
    % the periodic part of the Bloch function at each point in the Brillouin zone:
    %
    InnerSum = zeros(2*Qy, 2*Qx, Mylat, Mxlat);  % Fourier coefficients of Bloch function - preallocate
    OuterExp = zeros(2*Qy, 2*Qx, Mylat, Mxlat);  % Exponentials needed to compute Wannier functions - preallocate
    OuterExp2 = zeros(2*Qy, 2*Qx, Mylat, Mxlat);  % Exponentials needed to compute Wannier functions - preallocate
    counter_kxky = 0;       % Gets set to 1 after computing the first point inside the Brillouin zone.
                            % This is needed to ensure continuity of the eigenvectors of HM
                            % (which are Fourier coefficients of periodic part of the Bloch functions)
                            % with respect to {kx,ky}. See the calculation of "overlap" below.
                            %
    for iky = 1 : 2*Qy      % length(kygrid)

        k2 = nonrectang_y * kygrid(iky);    % wavenumbers (pseudo-momenta) of the Bloch functions
        eik2ylat = exp(1i *(kygrid(iky) * (ylat - siteY)) );
        eikylat = exp(1i *(kygrid(iky) * (ylat)) );

        for ikx = 1 : 2*Qx  % length(kxgrid) 

            k1 = nonrectang_x * kxgrid(ikx); 

            HM = zeros((2*Nx+1), (2*Ny+1), (2*Nx+1), (2*Ny+1));

            %  EXPLANATION OF SETUP OF THE MATRIX:
            %    Let iii, jjj be the Fourier indices of the Fourier transform G of the periodic part of
            %    the Bloch function along kx- and ky-directions, respectively. This G satisfies the equation:
            %    free_part(iii,jjj)*G(iii,jjj) + sum_{mmm,nnn} VF(iii-mmm, jjj-nnn)*G(mmm,nnn) = mu*G(iii,jjj),
            %    where:
            %    free_part(iii,jjj) = (k1 + iii*k0x)^2 + (k2 + jjj*k0y)^2 is the free part of the Hamiltonian;
            %    k1, k2 are the x-, y-wavenumbers in the Brillouin zone;
            %    k0{x,y} = 2pi/L_singlePeriod{x,y} = smallest wavenumbers of the periodic potential;
            %    VF(kkk,lll) are the Fourier coefficients of the potential;
            %    mu = energy (for the given (k1,k2)). 
            %    This equation can be rewritten as
            %    sum_{mmm,nnn} HM(mmm,nnn,iii,jjj) * G(mmm,nnn) = mu * G(iii,jjj),
            %    where:
            %    HM(mmm,nnn,iii,jjj) = 
            %      free_part(iii,jjj)*Kr_delta(iii,mmm)*Kr_delta(jjj,nnn) + VF(iii-mmm, jjj-nnn),
            %    and Kr_delta(a,b) is the Kronecker delta symbol. 
            %  
            %    The structure of HM is discussed further below, 
            %    where we flatten it into a  (2*Nx+1) x (2*Ny+1)  matrix. 
            % 
            for iii = 1 : (2*Nx+1)        % index along columns of HM (no relation to wavenumbers kx,ky)
                for jjj = 1 : (2*Ny+1)    % index along rows of z

                    for mmm = 1 : (2*Nx+1)
                        for nnn = 1 : (2*Ny+1)

                            if    iii - mmm + Nx+1 >= 1   &  iii - mmm + Nx+1 <= 2*Nx+1 ...
                               &  jjj - nnn + Ny+1 >= 1   &  jjj - nnn + Ny+1 <= 2*Ny+1

                                HM(mmm, nnn, iii, jjj) = - Vfourier(iii - mmm + Nx+1, jjj - nnn + Ny+1);  
                            end
                        end
                    end

                    % add the "free Hamiltonian" part to the diagonal
                    HM(iii, jjj, iii, jjj) = HM(iii, jjj, iii, jjj) + ...
                                             ( k1 + (iii - Nx - 1 )*norm_g1_nondim )^2 + ...
                                             ( k2 + (jjj - Ny - 1 )*norm_g2_nondim )^2 + ...
                                             2*( k1 + (iii - Nx - 1 )*norm_g1_nondim )* ...
                                               ( k2 + (jjj - Ny - 1 )*norm_g2_nondim )*cosine_g1g2; 

    %                 % This is for debugging only:
    %                 aux11(iii,jjj) = (k1 + (iii - N - 1 )*k0x)^2 + ...
    %                     (k2 + (jjj - N - 1 )*k0y)^2;
    %                 aux22(iii,jjj) = (k1 + (jjj - N - 1 )*k0x)^2 + ...
    %                     (k2 + (iii - N - 1 )*k0y)^2;
                end
            end

    %         % This is for debugging only:        
    %         [k1,k2]
    %         aux11-aux22
    %         pause

            % Flatten HM
            HM_2x2 = reshape(HM, [(2*Nx+1)*(2*Ny+1), (2*Nx+1)*(2*Ny+1)]) .';
            %%%%%%%%% HM_2x2 = reshape(permute(HM, [2,1,4,3]), [(2*N+1)^2, (2*N+1)^2]);
            % EXPLANATION: 
            % 1) The structure of HM can be viewed as an (2Nx+1) x (2Ny+1) array of (2Nx+1) x (2Ny+1) matrices.
            %    The (iii,jjj)-th matrix  HM(:,:,iii,jjj) multiplies the G(:,:) matrix *in the H'Adamard sense*
            %    (i.e., as sum(sum( HM(:,:,iii,jjj) .* G(:.:) )) ), and this yields the (iii,jjj)-th equation
            %    of the Fourier-transformed equation. 
            % 2) When Matlab flattens a matrix using  reshape,  it does so along the columns. 
            %    So, first, it would flatten each HM(:,:,iii,jjj) along its columns. 
            %    Then, it would put these columns in the order:
            %     (iii,jjj) = (1,1), (2, 1), ... (2Nx+1,1), (1,2), ... , (2Nx+1,2), (1,3), ... 
            % 3) We want each *column* of the flattened matrix HM to multiply an eigenvector.
            %    This is why we have the transpose after reshape. 
            %    If one does *not* use the transpose, the code will still work, as long as the potential V
            %    if real and symmetric, which would ensure that HM_2x2 is real and symmetric 
            %    (see the VERY IMPORTANT NOTE below). 
            %    However, just so that the meaning of each of the above steps would be transparent,
            %    we do use the transpose. 
            % 4) The rows of the so transposed matrix will be in the above order
            %    (1,1), (2, 1), ... (2Nx+1,1), (1,2), ... , (2Nx+1,2), (1,3), ... ,
            %    which will then be the order of the entries of the eigenvector Gsground, found below.
            %    Then, to obtain its (2Nx+1) x (2Ny+1) form, one simply uses  reshape,
            %    which will put it in the correctly ordered  (2Nx+1) x (2Ny+1)  matrix. 


            % 
            % VERY IMPORTANT NOTE: 
            % HM_2x2 = real & symmetric matrix !!!
            % This follows from the fact that the potential V in the stationary Schrodinger equation
            % is real and *even* in and in y; hence its Fourier coefficients are also real.
            % They are symmetric for the following reason (which is easiest to explain in 1D;
            % in 2D the situation is *probably* the same.
            % So, Fourier coefficients of V enter with indices V_{iii, mmm} = V_{iii - mmm}.
            % Then V_{mmm, iii} = V_{-(iii - mmm)}. Finally, the symmetry  V_jjj = V_{-jjj}
            % holds because V is real. 
            % Since HM is real and symmetric, its eigenvalues are real, and hence the eigenvectors 
            % computed by Matlab are also real. This will play role when we define their "overlap" below. 

    %         % For debugging only - view the structure and entries of HM_2x2:
    %         [k1  k2]
    %         figure(1077);
    %         mesh(HM_2x2)
    %         zlabel('HM_{2x2}')
    %         pause

            initial = ones((2*Nx+1)*(2*Ny+1), 1);
            options.v0 = initial;

            % FInd eigenvalues and eigenvectors of HM:
            [G, mu] = eigs(HM_2x2, Num_mu, 'smallestreal', options); %'smallestabs');
            [mus, indmu] = sort(diag(mu));
            Gs = G(:,indmu);
            Gsground = Gs(:,1);   % select the eigenvector corresponding to the lowest eigenvalue
                                  % (i.e., for the lowest energy band)
            %
            % Adjust the eigenvectors to depend continuously on the location inside the 
            % Brillouin zone.
            % EXPLANATION: 
            % Matlab finds real (see the VERY IMPORTANT NOTE above) eigenvectors up to the +/- sign.
            % We want this choice of sign to be continuous in {kx,ky}. Otherwise, when we integrate
            % the Bloch function (which is just the Fourier transform of the eigenvector)
            % over the Brillouin zone, we will get spurious (and wrong) results for the Wannier function.
            %
            if counter_kxky == 1  % Check that there has been a previous calculation of the eigenvector.
                overlap = Gsground.' * Gsground_old;
                if overlap < 0 
                    Gsground = - Gsground;
                end
            end
            %
            Gsground_grid = reshape(Gsground, [2*Nx+1, 2*Ny+1]);    % % % .'
            % This creates a grid Gij where i goes right along rows and j goes down along columns.

            counter_kxky = 1;      % Simply says that there has been a calculation for at least
                                   % one value of (k1,k2). This is needed for the "overlap"
                                   % computation above.
            Gsground_old = Gsground;


    %         % The block below is mostly for debugging and allows one to look at the eigenvectors
    %         % corresponding to *two* (as opposed to one) lowest eigenvalues.
    %         Gsground_1 = Gs(:,1);
    %         Gsground_2 = Gs(:,2);
    %         Gsground_1_grid = reshape(Gsground_1, [2*N+1, 2*N+1]);
    %         Gsground_2_grid = reshape(Gsground_2, [2*N+1, 2*N+1]);
    %         Gsground_grid = 1*Gsground_1_grid + 0*Gsground_2_grid;
    %
    %         % For debugging only (use EITHER figure 1001 or 1002):
    % %         [k1  k2]
    % %         mus
    %         figure(1001)
    %         subplot(2,2,1)
    %         mesh(Knvec/2, Kmvec/2, real(Gsground_1_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Re 1');
    %         subplot(2,2,2)
    %         mesh(Knvec/2, Kmvec/2, imag(Gsground_1_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Im 1');
    %         subplot(2,2,3)
    %         mesh(Knvec/2, Kmvec/2, real(Gsground_2_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Re 2');
    %         subplot(2,2,4)
    %         mesh(Knvec/2, Kmvec/2, imag(Gsground_2_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Im 2');
    %         %pause
    % %
    %         figure(1002)
    %         subplot(2,1,1)
    %         mesh(Knvec/2, Kmvec/2, real(Gsground_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Re Gs');
    %         subplot(2,1,2)
    %         mesh(Knvec/2, Kmvec/2, imag(Gsground_grid));
    %         xlabel('m'); ylabel('n'); zlabel('Im Gs');
    %         %pause


            % Save the eigenvalues of HM into a matrix which subsequently can be plotted
            % to show the band structure:
            Mu(iky,ikx,1:Num_mu) = mus(1:Num_mu);

            % InnerSum is the periodic part (in x and y) of the Bloch function for each (k1,k2).
            % OuterExp is used outside the loop over Brillouin zone to find the Wannier function.

            for ix = 1:Mxlat
                aux1 = Gsground_grid .* eiKmx(:,:, ix);
                aux2 = exp(1i * (kxgrid(ikx) * (xlat(ix) - siteX)) );
                aux3 = exp(1i * (kxgrid(ikx) * (xlat(ix))) );
               for iy = 1:Mylat
                   InnerSum(iky, ikx, iy, ix) = sum(sum( aux1 .* eiKny(:,:, iy) ));
                   OuterExp(iky, ikx, iy, ix) = aux2 * eik2ylat(iy); 
                   OuterExp2(iky, ikx, iy, ix) = aux3 * eikylat(iy); 
               end
            end

    %         % For debugging only - look at the periodic part of the Bloch function at one (k1,k2):
    %         figure(1003);
    %         % [k1,k2]
    %         InnerSum_aux(:,:)=InnerSum(iky, ikx, :, :);
    %         subplot(2,1,1);
    %         mesh(x,y,real(InnerSum_aux));
    %         xlabel('x'); ylabel('y'); zlabel('Real G');
    %         subplot(2,1,2);
    %         mesh(imag(InnerSum_aux));
    %         xlabel('x'); ylabel('y'); zlabel('Imag G');
    %         pause;

    %         % For debugging only - look (use TOP VIEW) at the exponentials at one (k1,k2):        
    %         figure(1004);
    %         [k1,k2]
    %         OuterExp_aux(:,:)=OuterExp(iky, ikx, :, :);
    %         subplot(2,1,1);
    %         mesh(x,y,real(OuterExp_aux)); view(2); axis('square')
    %         xlabel('x'); ylabel('y'); zlabel('Real Outer Exp');
    %         subplot(2,1,2);
    %         mesh(x,y,imag(OuterExp_aux)); view(2); axis('square')
    %         xlabel('x'); ylabel('y'); zlabel('Imag Outer Exp');
    %         pause;


            % loop over k1 and k2 ends here

        end
    end

    toc

    % % This is needed for debugging only - view the Bloch function as that of ky and y:
    % figure(1004)
    % aux(:,:) = abs(InnerSum(:,1,:,1));
    % mesh(aux);
    % 


    % ---
    % Compute Wannier function:
    disp('Computed the periodic part of the Bloch function and outer exponentials;')
    disp('Working on computing Wanner functions.')

    tic

    jacobian_a1a2_mat = [a1norm' a2norm'];  % [a1_nondim' a2_nondim']; %
    jacobian_a1a2_mat_inv = inv(jacobian_a1a2_mat);

    jacobian_a1a2 = det(jacobian_a1a2_mat);
    jacobian_g1g2 = det([g1norm' g2norm']);  % det([g1_nondim' g2_nondim']);

    %%%%%%%%%%%%%%%
    %%% - 10c --->
    jacobian_g1g2_dkxdky = jacobian_g1g2 * dkxgrid * dkygrid;

    Wannier = zeros(Mylat, Mxlat);
    for ix = 1:Mxlat
        for iy = 1:Mylat 
            Wannier(iy, ix) = sum(sum( OuterExp(:,:, iy, ix) .* ...
                                       InnerSum(:,:, iy, ix) )) * jacobian_g1g2_dkxdky;
        end
    end

    iiii = 12
    jjjj = 12

    BlochPsi = zeros(Mylat, Mxlat);
    for ix = 1:Mxlat
        for iy = 1:Mylat 
            BlochPsi(iy, ix) =  OuterExp2(iiii,jjjj, iy, ix) .* ...
                                       InnerSum(iiii,jjjj, iy, ix)  * jacobian_g1g2_dkxdky;
        end
    end

    Blochu = zeros(Mylat, Mxlat);
    for ix = 1:Mxlat
        for iy = 1:Mylat 
            Blochu(iy, ix) = InnerSum(iiii,jjjj, iy, ix);
        end
    end
    % Old normalization on rectangular grid
    % Wannier = Wannier * ((pi/Lxlat_singlePeriod)/Qx * (pi/Lylat_singlePeriod)/Qy) * ...
    %                     1/(4*pi);   % This latter normalization constant was worked out on paper.

    % Manual normalization of the wannier function
    %%%   Wannier_Norm = Wannier.*Wannier * (dxlat_singlePeriod*dylat_singlePeriod);

    C = sum(sum(Wannier.^2)) * jacobian_a1a2 * (dxlat_singlePeriod*dylat_singlePeriod);
    %%% <--- 10c -
    %%%%%%%%%%%%%%%

    Wannier = Wannier/sqrt(C);

    % ---
    % Comapare our 2D results with Adrian's 1D result:
    %
    % Loads 1d Wannier function:
    % Wannier_1D = readmatrix('wannier10.txt');   % This doesn't work with Matlab 2017b. 
    % Wannier_1D = dlmread('wannier10.txt');
    % C_1D = 1/sqrt(dx_singlePeriod * sum(Wannier_1D.^2));  % normalization constant by brute force
    % Wannier_1D_Norm = C_1D * Wannier_1D;
    % [W1x, W1y] = meshgrid(Wannier_1D_Norm, Wannier_1D_Norm);
    % Wannier_combined = W1x.*W1y;

    % % For debugging only:
    % %disp('normalization constant');
    % C_2D = 1/sqrt(dx_singlePeriod * dy_singlePeriod * sum(sum(abs(Wannier).^2)));
    % Wannier_2D_Norm = C_2D * Wannier;
    Wannier_2D_Norm = Wannier;

    pimc_dens_lat_arr = readmatrix("pimc_dens_lat.txt");
    pimc_psi_lat_arr = sqrt(pimc_dens_lat_arr);

    C_w = sum(sum(pimc_psi_lat_arr.^2)) * jacobian_a1a2 * (dxlat_singlePeriod*dylat_singlePeriod);

    pimc_wannier = pimc_psi_lat_arr/sqrt(C_w);

    % Only for using PIMC dens as Wannier Function
    %Wannier_2D_Norm = reshape(abs(pimc_wannier), [288, 288]);

    toc

    V_graphene_lat = V_multiplier*dlmread(strcat('V_graphene_lat', ...
                                                  strain_strings(iss) ,'.txt'));
    V_graphene_lat_mod = V_graphene_lat(1:end-1, 1:end-1);

    V_potential_graphene = sum(sum(V_graphene_lat_mod.*Wannier_2D_Norm.^2))*...
                           jacobian_a1a2 * dylat_singlePeriod * dxlat_singlePeriod;


    % ----------------------------------------------------------------------------------
    % Integral for the hopping coefficent using tight binding approx,
    %%%%%%%%%%%%%%%
    %%% - 10c ---> 
    % computed via the integration of the band energy over the BZ:
    %%% <--- 10c -
    %%%%%%%%%%%%%%%

    t2_sum_a1 = zeros(2*Qy, 2*Qx);
    t2_sum_a2 = zeros(2*Qy, 2*Qx);

    %%%%%%%%%%%%%%%
    %%% - 10c --->
    % jacobian_g1g2_dkxdky has been defined above 
    % jacobian_g1g2_dkdk = jacobian_g1g2 * (k0x/2)/Qx * (k0y/2)/Qy;

    for iky = 1 : 2*Qy      % length(kygrid)
        for ikx = 1 : 2*Qx
            t2_sum_a1(iky, ikx) = Mu(iky,ikx,1) * exp(-1i *...
                                 (Lxlat_singlePeriod * kxgrid(ikx) + ...  
                                  ... % previously used norm(a1) instead of Lxlat_singlePeriod 
                                 0*kygrid(iky)));
            t2_sum_a2(iky, ikx) = Mu(iky,ikx,1) * exp(-1i *...
                                 (0*kxgrid(ikx) + ...
                                 Lylat_singlePeriod * kygrid(iky)));
                                 ... % previously used norm(a2) instead of Lylat_singlePeriod 
        end
    end
    %%% <--- 10c -
    %%%%%%%%%%%%%%%

    t2_a1 = 1/A_BZ * sum(sum(t2_sum_a1)) * jacobian_g1g2_dkxdky * nonrectang_x*nonrectang_y;
    t2_a2 = 1/A_BZ * sum(sum(t2_sum_a2)) * jacobian_g1g2_dkxdky * nonrectang_x*nonrectang_y;
    % The "1/A_BZ" term is required by the analytical formula.
    % The "nonrectang_x*nonrectang_y" term is required because our dkx and dky are "shorter" 
    % than the reciprocal vectors by the factors nonrectang_x and nonrectang_y. 

    disp('t computed via integration over BZ:')
    [ t2_a1 
      t2_a2 ]


    % ----------------------------------------------------------------------------------



    tic


    for iy = 1:Mylat
        Dx_Wannier_2D_Norm(iy,:) = [diff(Wannier_2D_Norm(iy,:)) 0]./dxlat_singlePeriod;
    end

    for ix = 1:Mylat
        Dy_Wannier_2D_Norm(:,ix) = [diff(Wannier_2D_Norm(:, ix)); 0]./dylat_singlePeriod;
    end

    for iy = 1:Mylat
        D2x_Wannier_2D_Norm(iy,:) = [diff(Dx_Wannier_2D_Norm(iy,:)) 0]./dxlat_singlePeriod;
    end

    for iy = 1:Mylat
        DxDy_Wannier_2D_Norm(iy,:) = [diff(Dy_Wannier_2D_Norm(iy,:)) 0]./dxlat_singlePeriod;
    end

    for ix = 1:Mylat
        DyDx_Wannier_2D_Norm(:,ix) = [diff(Dx_Wannier_2D_Norm(:, ix)); 0]./dylat_singlePeriod;
    end

    for ix = 1:Mylat
        D2y_Wannier_2D_Norm(:,ix) = [diff(Dy_Wannier_2D_Norm(:, ix)); 0]./dylat_singlePeriod;
    end

    %%%%%%%%%%%%%
    % Alternative method to compute derivatives is via Discrete Fourier Transform:
    %
    dwavenumber_x = 2*pi/Lxlat; 
    dwavenumber_y = 2*pi/Lylat;
    wavenumber_x = dwavenumber_x* [0 : Mxlat/2-1  -Mxlat/2 : -1];
    wavenumber_y = dwavenumber_y* [0 : Mylat/2-1  -Mylat/2 : -1];
    wavenumber_x_sq = wavenumber_x.^2;
    wavenumber_y_sq = wavenumber_y.^2;
    %

    D2x_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));
    Dx_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));
    D2y_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));
    Dy_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));
    DxDy_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));
    DyDx_ALT_Wannier_2D_Norm = zeros(size(Wannier_2D_Norm));

    for iy = 1:Mylat
        fftx_Wannier = fft(Wannier_2D_Norm(iy,:));
        D2x_ALT_Wannier_2D_Norm(iy,:) = ifft( -wavenumber_x_sq .* fftx_Wannier );
        Dx_ALT_Wannier_2D_Norm(iy,:)  = 1i*ifft( wavenumber_x .* fftx_Wannier );
    end
    for ix = 1:Mxlat
        ffty_Wannier = fft(Wannier_2D_Norm(:,ix));
        D2y_ALT_Wannier_2D_Norm(:,ix) = ifft( -wavenumber_y_sq' .*  ffty_Wannier);
        Dy_ALT_Wannier_2D_Norm(:,ix)  = 1i*ifft( wavenumber_y' .* ffty_Wannier );
        DxDy_ALT_Wannier_2D_Norm(:,ix) = 1i*ifft( wavenumber_y' .* fft(Dx_ALT_Wannier_2D_Norm(:,ix)) );
    end
    for iy = 1:Mylat
        fftx_Wannier = fft(Wannier_2D_Norm(iy,:));
        DyDx_ALT_Wannier_2D_Norm(iy,:) = 1i*ifft( wavenumber_x .* fft(Dy_ALT_Wannier_2D_Norm(iy,:)) );
    end


    % Translate the position of the Wannier function
    Wannier_a1 = zeros(size(Wannier_2D_Norm));
    %     D2x_ALT_Wannier_a1_Norm = zeros(size(Wannier_2D_Norm));
    %     Dx_ALT_Wannier_a1_Norm = zeros(size(Wannier_2D_Norm));
    %     D2y_ALT_Wannier_a1_Norm = zeros(size(Wannier_2D_Norm));
    %     Dy_ALT_Wannier_a1_Norm = zeros(size(Wannier_2D_Norm));
    %     DxDy_ALT_Wannier_a1_Norm = zeros(size(Wannier_2D_Norm));

    for ix = 1 + Mxlat_singlePeriod : Mxlat
        for iy = 1 : Mylat 
            Wannier_a1(iy,ix) = Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
    %             D2x_ALT_Wannier_a1_Norm(iy,ix) = D2x_Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
    %             Dx_ALT_Wannier_a1_Norm(iy,ix) = Dx_Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
    %             D2y_ALT_Wannier_a1_Norm(iy,ix) = D2y_Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
    %             Dy_ALT_Wannier_a1_Norm(iy,ix) = Dy_Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
    %             DxDy_ALT_Wannier_a1_Norm(iy,ix) = DxDy_Wannier_2D_Norm(iy,ix - Mxlat_singlePeriod);
        end
    end

    % 2nd shifted wannier functions
    Wannier_a2 = zeros(size(Wannier_2D_Norm));
    Wannier_a7 = zeros(size(Wannier_2D_Norm));

    for ix = 1 : Mxlat
        for iy = 1 + Mylat_singlePeriod : Mylat 
            Wannier_a2(iy,ix) = Wannier_2D_Norm(iy - Mylat_singlePeriod, ix);
            Wannier_a7(iy,ix) = Wannier_a1(iy - Mylat_singlePeriod, ix);
        end
    end
    % 
    % 3rd and 4th shifted wannier functions
    Wannier_a6 = zeros(size(Wannier_2D_Norm));
    Wannier_a5 = zeros(size(Wannier_2D_Norm));

    for ix = 1 : Mxlat
        for iy = 1 : Mylat - Mylat_singlePeriod
            Wannier_a6(iy,ix) = Wannier_a1(iy + Mylat_singlePeriod, ix);
            Wannier_a5(iy,ix) = Wannier_2D_Norm(iy + Mylat_singlePeriod, ix);
        end
    end

    Wannier_a3 = zeros(size(Wannier_2D_Norm));
    Wannier_a4 = zeros(size(Wannier_2D_Norm));

    for ix = 1 : Mxlat - Mxlat_singlePeriod
        for iy = 1 : Mylat
            Wannier_a3(iy,ix) = Wannier_a2(iy, ix + Mxlat_singlePeriod);
            Wannier_a4(iy,ix) = Wannier_2D_Norm(iy,ix + Mxlat_singlePeriod);
        end
    end
    %    
    %     % Translate the position of the Wannier function
    %     Wannier_a7 = zeros(size(Wannier_2D_Norm));
    %     Wannier_a12 = zeros(size(Wannier_2D_Norm));
    % 
    %     for ix = 1 + Mxlat_singlePeriod : Mxlat
    %         for iy = 1 : Mylat 
    %             Wannier_a7(iy,ix) = Wannier_a2(iy,ix - Mxlat_singlePeriod);
    %             Wannier_a12(iy,ix) = Wannier_a6(iy,ix - Mxlat_singlePeriod);
    %         end
    %     end
    %     
    %     Wannier_a8 = zeros(size(Wannier_2D_Norm));
    % 
    %     for ix = 1 : Mxlat
    %         for iy = 1 + Mylat_singlePeriod : Mylat 
    %             Wannier_a8(iy,ix) = Wannier_a3(iy - Mylat_singlePeriod, ix);
    %         end
    %     end
    % 
    %     Wannier_a11 = zeros(size(Wannier_2D_Norm));
    %     Wannier_a10 = zeros(size(Wannier_2D_Norm));
    % 
    %     for ix = 1 : Mxlat
    %         for iy = 1 : Mylat - Mylat_singlePeriod
    %             Wannier_a11(iy,ix) = Wannier_a6(iy + Mylat_singlePeriod, ix);
    %             Wannier_a10(iy,ix) = Wannier_a4(iy + Mylat_singlePeriod, ix);
    %         end
    %     end
    % 
    %     Wannier_a9 = zeros(size(Wannier_2D_Norm));
    % 
    %     for ix = 1 : Mxlat - Mxlat_singlePeriod
    %         for iy = 1 : Mylat
    %             Wannier_a9(iy,ix) = Wannier_a3(iy, ix + Mxlat_singlePeriod);
    %         end
    %     end

    % --------------------------------------------------------------------------
    a1a2jacobian_dxdy_sq = jacobian_a1a2^2 * dylat_singlePeriod^2 * dxlat_singlePeriod^2;

    He_He = dlmread('He_Aziz95.dat', '', 0, 0);
    r_He_He = [space_normalization.*He_He(:,1); ...
               1.00001*space_normalization.*He_He(end,1) ;1000];
    V_He_He = energy_normalization*[He_He(:,2); 0; 0];

    J_exp = 0 .* r_He_He;

    for ic = 1:length(J_exp)
        if r_He_He(ic) < cutoff
            J_exp(ic) = 1;
        else
            J_exp(ic) = 1;
        end
    end

    reducesize_forU = 2;

    Mxlat_rsforU = 2*round( Mxlat / (2*reducesize_forU) );
    Mylat_rsforU = 2*round( Mylat / (2*reducesize_forU) );

    Lxlat_rsforU = Lxlat/reducesize_forU;
    Lylat_rsforU = Lylat/reducesize_forU;


    % creates array of distances from point Z_He_He above x-y plane to each
    % lattice point
    rlat_arr = zeros(2*Mylat_rsforU-1, 2*Mxlat_rsforU-1);  % ---> Jiang, it was zeros(Mylat, Mxlat); 
                                                            %      did you mean zeros(2*Mylat-1, 2*Mxlat-1);

    for iy = -(Mylat_rsforU) : (Mylat_rsforU - 1)
        for ix = -(Mxlat_rsforU) : (Mxlat_rsforU - 1)
            rlat_arr(iy + Mylat_rsforU + 1, ix + Mxlat_rsforU + 1) = sqrt(...
                (ix*a1(1)/Mxlat_singlePeriod + iy * a2(1)/Mylat_singlePeriod)^2 + ...
                (ix*a1(2)/Mxlat_singlePeriod + iy * a2(2)/Mylat_singlePeriod)^2);

            xlat_rsforU(iy + Mylat_rsforU + 1, ix + Mxlat_rsforU + 1) = ...
                ix * a1(1) + iy * a2(1);
            ylat_rsforU(iy + Mylat_rsforU + 1, ix + Mxlat_rsforU + 1) = ...
                ix * a1(2) + iy * a2(2);
        end
    end


    rlat_vec = reshape(rlat_arr, [1 ,(2*Mylat_rsforU)*(2*Mxlat_rsforU)]);

    % Calculates the V values to be accessed
    V_He_He_vec = interp1(r_He_He, V_He_He, rlat_vec);
    J_exp_vec = interp1(r_He_He, J_exp, rlat_vec);

    V_He_He_arr = zeros(2*Mylat_rsforU);
    J_exp_arr = zeros(2*Mylat_rsforU);

    V_He_He_arr = reshape(V_He_He_vec, [2*Mylat_rsforU, 2*Mxlat_rsforU]);
    J_exp_arr = reshape(J_exp_vec, [2*Mylat_rsforU, 2*Mxlat_rsforU]);

    V_Cutoff_arr = V_He_He_arr.*J_exp_arr;

    Wannier_0_sq = Wannier_2D_Norm.^2;
    Wannier_a1_sq = Wannier_a1.^2;
    Wannier_a2_sq = Wannier_a2.^2;
    Wannier_a3_sq = Wannier_a3.^2;
    Wannier_a7_sq = Wannier_a7.^2;

    fftwf0_sq = fft2(Wannier_0_sq);
    fftwf1_sq = fft2(Wannier_a1_sq);
    fftwf2_sq = fft2(Wannier_a2_sq);
    fftwf3_sq = fft2(Wannier_a3_sq);
    fftwf7_sq = fft2(Wannier_a7_sq);
    fftVHeHe = fft2(V_Cutoff_arr);
    % for testing
    fftones = fft2(ones(2*Mylat_rsforU));

    V_BH_ALT_a1 = sum(sum( Wannier_0_sq .* fftshift( ifft2( fftwf1_sq .* fftVHeHe ) ) ) ) * a1a2jacobian_dxdy_sq;
    V_BH_ALT_a2 = sum(sum( Wannier_0_sq .* fftshift( ifft2( fftwf2_sq .* fftVHeHe ) ) ) ) * a1a2jacobian_dxdy_sq;
    V_BH_ALT_a3 = sum(sum( Wannier_0_sq .* fftshift( ifft2( fftwf3_sq .* fftVHeHe ) ) ) ) * a1a2jacobian_dxdy_sq;
    V_BH_ALT_a7 = sum(sum( Wannier_0_sq .* fftshift( ifft2( fftwf7_sq .* fftVHeHe ) ) ) ) * a1a2jacobian_dxdy_sq;
    % for testing
    Norm_test_conv = sum(sum( Wannier_0_sq .* fftshift( ifft2( fftwf0_sq .* fftones ) ) ) ) * a1a2jacobian_dxdy_sq;

    % Integral for the hopping coefficent using wannier functions
    V_graphene_lat = V_multiplier*dlmread(strcat('V_graphene_lat', ...
                                                  strain_strings(iss) ,'.txt'));

    t_sum_a1 = zeros(Mylat, Mxlat);
    t_sum_a2 = zeros(Mylat, Mxlat);
    t_sum_a3 = zeros(Mylat, Mxlat);


    a1a2jacobian_dxdy = jacobian_a1a2 * dylat_singlePeriod * dxlat_singlePeriod;

    % Coefficents of the Lapacian
    %%%%%%%%%%%%%%%
    %%% - 10c --->
    c_laplace_xlxl = - (jacobian_a1a2_mat_inv(1,1)^2 + jacobian_a1a2_mat_inv(1,2)^2);
    c_laplace_ylyl = - (jacobian_a1a2_mat_inv(2,1)^2 + jacobian_a1a2_mat_inv(2,2)^2);
    c_laplace_xlyl = - 2*(jacobian_a1a2_mat_inv(1,1)*jacobian_a1a2_mat_inv(2,1) + ...
                          jacobian_a1a2_mat_inv(1,2)*jacobian_a1a2_mat_inv(2,2));
    % c_laplace_xlxl = -norm(a1_nondim)^2*(jacobian_a1a2_mat_inv(1,1)^2 + jacobian_a1a2_mat_inv(1,2)^2);
    % c_laplace_ylyl = -norm(a1_nondim)^2*(jacobian_a1a2_mat_inv(2,1)^2 + jacobian_a1a2_mat_inv(2,2)^2);
    % c_laplace_xlyl = -2*norm(a1_nondim)^2*(jacobian_a1a2_mat_inv(1,1)*jacobian_a1a2_mat_inv(2,1) + ...
    %                  jacobian_a1a2_mat_inv(1,2)*jacobian_a1a2_mat_inv(2,2));
    %%% <--- 10c -
    %%%%%%%%%%%%%%%

    for iy = 1:Mylat
        for ix = 1:Mxlat
            %%%%%%%%%%%%%%%
            %%% - 10c --->
            aux1 = ( (c_laplace_xlxl * D2x_Wannier_2D_Norm(iy, ix) + ...   % was a wrong parenthesis
                      c_laplace_xlyl * DxDy_Wannier_2D_Norm(iy, ix) + ...
                      c_laplace_ylyl * D2y_Wannier_2D_Norm(iy, ix) ) + ...
                      V_graphene_lat(iy, ix) * ...
                      Wannier_2D_Norm(iy, ix) ) * ...
                   a1a2jacobian_dxdy;
            aux2 = ( (c_laplace_xlxl * D2x_Wannier_2D_Norm(iy, ix) + ...   % parenthesis was okay here
                      c_laplace_xlyl * DyDx_Wannier_2D_Norm(iy, ix) + ...
                      c_laplace_ylyl * D2y_Wannier_2D_Norm(iy, ix)) + ...
                      V_graphene_lat(iy, ix) * ...
                      Wannier_2D_Norm(iy, ix) ) * ...
                   a1a2jacobian_dxdy;

            aux1_ALT = ( (c_laplace_xlxl * D2x_ALT_Wannier_2D_Norm(iy, ix) + ...   % was a wrong parenthesis
                      c_laplace_xlyl * DxDy_ALT_Wannier_2D_Norm(iy, ix) + ...
                      c_laplace_ylyl * D2y_ALT_Wannier_2D_Norm(iy, ix) ) + ...
                      V_graphene_lat(iy, ix) * ...
                      Wannier_2D_Norm(iy, ix) ) * ...
                      a1a2jacobian_dxdy;

            auxarr(iy, ix) = aux1_ALT;

            t_sum_a1(iy, ix) = Wannier_a1(iy, ix) .*...    % used "*" instead of the correct ".*"
                               aux1;
            t_sum_a2(iy, ix) = Wannier_a2(iy, ix) .*...    % used "*" instead of the correct ".*"
                               aux2;
            t_sum_a3(iy, ix) = Wannier_a3(iy, ix) .*...    % used "*" instead of the correct ".*"
                               aux2;

            t_ALT_sum_a1(iy, ix) = Wannier_a1(iy, ix) .* aux1_ALT;
            t_ALT_sum_a2(iy, ix) = Wannier_a2(iy, ix) .* aux1_ALT;
            t_ALT_sum_a3(iy, ix) = Wannier_a3(iy, ix) .* aux1_ALT;

            forbenchmarking1_sum(iy, ix) = Wannier_2D_Norm(iy, ix) * Wannier_2D_Norm(iy, ix);
            forbenchmarking2_sum(iy, ix) = Wannier_2D_Norm(iy, ix) * aux1;
            %%% <--- 10c -
            %%%%%%%%%%%%%%%

        end
    end

    t_a1 = sum(sum(t_sum_a1));
    t_a2 = sum(sum(t_sum_a2));
    t_a3 = sum(sum(t_sum_a3));
    t_ALT_a1 = sum(sum(t_ALT_sum_a1));
    t_ALT_a2 = sum(sum(t_ALT_sum_a2));
    t_ALT_a3 = sum(sum(t_ALT_sum_a3));
    forbenchmarking1 = sum(sum( forbenchmarking1_sum ))*a1a2jacobian_dxdy;
    forbenchmarking2 = sum(sum( forbenchmarking2_sum ));

    %     V_BH_a1_array(iss) = V_BH_ALT_a1/energy_normalization;
    %     V_BH_a2_array(iss) = V_BH_ALT_a2/energy_normalization;
    %     V_BH_a3_array(iss) = V_BH_ALT_a3/energy_normalization;
    %     V_BH_a7_array(iss) = V_BH_ALT_a7/energy_normalization;
    t_a1_array(iss) = t_a1;
    t_a2_array(iss) = t_a2;
    t_a3_array(iss) = t_a3;
    t_a1_ALT_array(iss) = t_ALT_a1;
    t_a2_ALT_array(iss) = t_ALT_a2;
    t_a3_ALT_array(iss) = t_ALT_a3;

    disp('t computed via integration over space:')
    [ t_a1
      t_a2 
      t_a3 
      t_ALT_a1
      t_ALT_a2
      t_ALT_a3
    %       V_BH_ALT_a1
    %       V_BH_ALT_a2
    %       V_BH_ALT_a3
    %       V_BH_ALT_a7
      forbenchmarking1 
      forbenchmarking2 ]
    % --------------------------------------------------------------------------
    
    toc 
end
 
% realWF = readmatrix("lineardensity-T-reduce-t-0.00313-u-+450.000-L-005.500-6eda2407-b71d-4de0-9f71-3bf32f629ffc.dat")
% 
% z_real = realWF(:,1) + 2.75
% rho_real = realWF(:,2)
% 
% rho_real = rho_real/(sum(rho_real)*(z_real(2)-z_real(1)))
% 
% rho_zs = interp1(z_real, rho_real, zs)
% 
% deltaV_eff = trapz(deltaV_arr.*rho_zs)*0.1
% 
% t_eff = trapz(t_a1_array.*rho_zs)*0.1

% x_cart_rec = zeros(1, 71);
% 
% for ix = 0:70
%     x_cart_rec(ix+1) = ((3+ix)*a1(1) + (73-ix)*a2(1))/space_normalization/24;
% end
% 
% y_cart_rec = zeros(1, 71);
% 
% for iy = 0:70
%     y_cart_rec(iy+1) = ((iy-70)*a1(2) + (iy)*a2(2))/space_normalization/24;
% end
% 
% Wannier_sq_rec = zeros(71, 71);
% Wannier_12pd_red = zeros(71, 71);
% Wannier_sq_arr = zeros(1, 71*71);
% Wannier_pd_arr = zeros(1, 71*71);
% x_arr = zeros(1, 71*71);
% y_arr = zeros(1, 71*71);
% 
% ia = 1
% for ix = 1:71
%     for iy = 1:71
%         sy = 73 + iy - ix;
%         sx = ix + iy + 1;
%         Wannier_sq_rec(iy, ix) = Wannier_rsforU(sy, sx).^2*space_normalization^2; 
%         Wannier_12pd_red(iy, ix) = prob_dens_r2(sy, sx)*space_normalization^2; 
%     end
% end
% 
% for ix = 1:71
%     Wannier_sq_rec(:, ix) = interp1(y_cart_rec, Wannier_sq_rec(:, ix), x_cart_rec); 
%     Wannier_12pd_red(:, ix) = interp1(y_cart_rec, Wannier_12pd_red(:, ix), x_cart_rec); 
% end
% 
% y_cart_rec = x_cart_rec;
% 
% for ix = 1:71
%     for iy = 1:71
%         Wannier_sq_arr(ia) = Wannier_sq_rec(iy, ix); 
%         Wannier_pd_arr(ia) = Wannier_12pd_red(iy, ix); 
%         x_arr(ia) = x_cart_rec(ix);
%         y_arr(ia) = y_cart_rec(iy);
%         ia = ia + 1;
%     end
% end
% 
% [x_rec_arr, y_rec_arr] = meshgrid(x_cart_rec, y_cart_rec);
% 
% out1 = [x_arr' y_arr' Wannier_sq_arr' Wannier_pd_arr'];
% 
% writematrix(real(out1), 'cartesian_wannier.txt')
% writematrix(real(Wannier_sq_rec), 'wannier_squared.txt')
% writematrix(real(Wannier_12pd_red), 'jastrow_sp_density.txt')
% 
% coord = [x_cart_rec' y_cart_rec'];
% writematrix(real(x_cart_rec), 'xy_coords.txt')
%  
% % CHECK: DOES JUPYTER PRODUCE RIGHT POTENTIAL???
% % Loading the potential from file
% 
V_graphene = dlmread('V_graphenea00.txt');
Julia_grid = dlmread('Julia_grida00.txt');
x_julia = Julia_grid(1,:);
y_julia = Julia_grid(2,:);
% V_sign = -1;
% V = V_sign .* V_graphene;
% V=-15./(1+16*cos(X_singlePeriod).^2.*cos(Y_singlePeriod).^2); % Jianke's potential
% V = -10.*(sin(X_singlePeriod).^2 + sin(Y_singlePeriod).^2 - 0);
% V = 1*sin(X_singlePeriod).^2;

% ---
% Plot results:

% Plot the potential:

% Real stands for real space

figure(110);
mesh(x_julia, y_julia, V_graphene);
xlabel('x'); ylabel('y'); zlabel('V')

% Plot energy bands:
g1g2matr = [g1'/norm(g1) g2'/norm(g2)];

[kxgrid_mesh, kygrid_mesh] = meshgrid(kxgrid, kygrid);
kxgird_vec = reshape(kxgrid_mesh, [1 numel(kxgrid_mesh)]);
kygird_vec = reshape(kygrid_mesh, [1 numel(kygrid_mesh)]);

kxkymatr = g1g2matr*[kxgird_vec; kygird_vec];

tReciprocal = delaunay(kxkymatr(1,:)', kxkymatr(2,:)');
% 
figure(111);
for imu = 1:Num_mu
%     mesh(kxgrid, kygrid, Mu(:,:,imu));
    trisurf(tReciprocal, kxkymatr(1,:)', kxkymatr(2,:)', ...
            reshape(Mu(:,:,imu), [1 numel(Mu(:,:,imu))]));
    hold on
end
axis([-2, 2, -2, 2, -18, -10]);
axis('square');
set(gca,'FontSize',12)
xlabel('x');
ylabel('y');
hold off;

% subplot(2,1,1)
% mesh(kxgrid, kygrid, Mu(:,:,1));
% xlabel('kx'); ylabel('ky'); zlabel('\mu_1');
% subplot(2,1,2)
% mesh(kxgrid, kygrid, Mu(:,:,2));
% xlabel('kx'); ylabel('ky'); zlabel('\mu_2');

% Plot the abs value of Wannier function:
% 
a1a2matr = [a1norm' a2norm'];
 
[xlat_mesh, ylat_mesh] = meshgrid(xlat, ylat);
xlat_vec = reshape(xlat_mesh, [1 numel(xlat_mesh)]);
ylat_vec = reshape(ylat_mesh, [1 numel(ylat_mesh)]);

xymatr = a1a2matr*[xlat_vec; ylat_vec];

tReal = delaunay(xymatr(1,:)'/space_normalization, xymatr(2,:)'/space_normalization);

figure(112);
trisurf(tReal, xymatr(1,:)', xymatr(2,:)', ...
        reshape(abs(Wannier_2D_Norm), [1 numel(abs(Wannier_2D_Norm))]),...
        'edgecolor','none');
hold on
trisurf(tReal, xymatr(1,:)', xymatr(2,:)', ...
        reshape(abs(Wannier_a1), [1 numel(abs(Wannier_a1))]),...
        'edgecolor','none');
hold on
trisurf(tReal, xymatr(1,:)', xymatr(2,:)', ...
        reshape(abs(Wannier_a3), [1 numel(abs(Wannier_a3))]),...
        'edgecolor','none');
hold off
xlabel('x'); ylabel('y'); zlabel('|Wannier| renormed')
 axis([-10, 10, -10, 10, 0, 0.8]);
% axis('square')
% mesh(xlat/Lxlat_singlePeriod, ylat/Lylat_singlePeriod, abs(Wannier_2D_Norm));
% hold on
% mesh(xlat/Lxlat_singlePeriod, ylat/Lylat_singlePeriod, abs(Wannier_a1));
% hold on
% mesh(xlat/Lxlat_singlePeriod, ylat/Lylat_singlePeriod, abs(Wannier_a2));
% hold on
% axis([-40, 40, -40, 40, 0, 0.5]);
% axis('square')

figure(121);
trisurf(tReal, xymatr(1,:)', xymatr(2,:)', ...
        reshape(abs(auxarr), [1 numel(abs(auxarr))]),...
        'edgecolor','none');
xlabel('x'); ylabel('y'); zlabel('|Wannier| renormed')
axis([-10, 10, -10, 10, 0, 0.2]);

% figure(3)
% mesh(xlat_rsforU(1:144,1:144), ylat_rsforU(1:144,1:144), abs(Wannier_rsforU))
% hold on 
% mesh(xlat_rsforU(1:144,1:144), ylat_rsforU(1:144,1:144), abs(Wannier_rsforU_a1))
% hold on
% mesh(xlat_rsforU(1:144,1:144), ylat_rsforU(1:144,1:144), abs(Wannier_rsforU_a3))
% hold off
% axis([-200,200,-600,-200,0,1])


% Plot the abs value of Wannier function:
% figure(122);
% %mesh(xlat/Lxlat_singlePeriod, ylat/Lylat_singlePeriod, log10(abs(Wannier_2D_Norm)));
% trisurf(tReal, xymatr(1,:)', xymatr(2,:)', ...
%         reshape(log10(abs(Wannier_2D_Norm)), [1 numel(abs(Wannier_2D_Norm))]),...
%         'edgecolor','none');
% %hold on
% xlabel('x'); ylabel('y'); zlabel('log |Wannier| renormed')
% % axis([-40, 40, -40, 40, -12, 0]);
% colorbar('eastoutside');
% axis('equal')


% Make a Gamma-M-K plot.
%
% 1) define a mesh that is twice as fine in kx:
finefactor = 2;   % <--- THIS NUMBER MUST NOT BE CHANGED!!!!!!!!!!!!
Qxfiner = Qx*finefactor;
kxgrid_finer = ( -(k0x/2) : (k0x/2)/Qxfiner : (k0x/2) - (k0x/2)/Qxfiner )  + (k0x/2)/Qxfiner; 
[kxgrid_finerx_mesh, kygrid_finerx_mesh] = meshgrid(kxgrid_finer, kygrid);
% 2) interpolate Mu on this finer mesh:
for imu = 1:Num_mu
    Mu_finerx(:,:,imu) = interp2(kxgrid_mesh, kygrid_mesh, Mu(:,:,imu), ...
                                 kxgrid_finerx_mesh, kygrid_finerx_mesh);
end
% 3) create segments of the path
%   a) Gamma -> M 
for imu = 1:Num_mu
    Mu_pathGM(:, imu) = Mu_finerx(Qy, 2*Qx : 4*Qx-1, imu);
end
% Gets proper lengths along the path
L_pathGM = sqrt(3)/2*[0:2*Qx-1;].*dxlat_singlePeriod;
%   b) M -> K 
for imu = 1:Num_mu
    for iii = 0 : 2*Qx/3-1
        Mu_pathMK(iii+1, imu) = Mu_finerx(Qy - iii, 4*Qx - iii, imu);
    end
end
% Gets proper lengths along the path
L_pathMK(1) = sqrt(3)/2*2*Qx*dxlat_singlePeriod;
for iii = 1 : 2*Qx/3-1
    L_pathMK(iii+1) = L_pathMK(iii) + ...
                      sqrt(dxlat_singlePeriod^2 + dylat_singlePeriod^2); 
end
%   c) K -> G
for imu = 1:Num_mu
    for iii = 0 : 2*Qx/3
        Mu_pathKG(iii+1, imu) = Mu_finerx(Qy/3 + iii, 10*Qx/3 - 2*iii, imu);
    end
end
L_pathKG(1) = L_pathMK(end) + sqrt(dxlat_singlePeriod^2 + dylat_singlePeriod^2); 
% Gets proper lengths along the path
for iii = 1 : 2*Qx/3
    L_pathKG(iii+1) = L_pathKG(iii) + ...
                      2*sqrt(4*dxlat_singlePeriod^2 + dylat_singlePeriod^2)/sqrt(3);
end

figure(113);
for imu = 1:Num_mu
    plot( [L_pathGM  L_pathMK  L_pathKG], ...
         ([Mu_pathGM(:, imu);  Mu_pathMK(:, imu);  Mu_pathKG(:, imu)]));
    hold on
end
% xline(L_pathMK(1));
% xline(L_pathKG(1));
lowerlim_GMK = min( [Mu_pathGM(:, 1);  Mu_pathMK(:, 1);  Mu_pathKG(:, 1)] );
upperlim_GMK = max( [Mu_pathGM(:, Num_mu);  Mu_pathMK(:, Num_mu);  Mu_pathKG(:, Num_mu)] );
line( [L_pathMK(1) L_pathMK(1)],  [lowerlim_GMK  upperlim_GMK] );
line( [L_pathKG(1) L_pathKG(1)],  [lowerlim_GMK  upperlim_GMK] );
xticks([0 L_pathMK(1) L_pathKG(1) L_pathKG(end)]);
xticklabels({'\Gamma', 'M', 'K', '\Gamma'})
hold off
ylim([lowerlim_GMK-0.01,  upperlim_GMK+0.01])
xlim([L_pathGM(1), L_pathKG(end)])
ylabel('E / E_R')
            
t_fromGMKspaghetti = (Mu_pathKG(1,1) - Mu_pathKG(end,1))/9;

disp("dV")
disp((max(max(V_graphene_lat)) - min(min(V_graphene_lat)))/energy_normalization)

Mus = [Mu_pathGM(1, imu);  Mu_pathMK(1, imu);  Mu_pathKG(1, imu)];

%disp("1st Bandwidth")
%disp((max(Mu_pathGM(:, 1))-min(Mu_pathGM(:, 1)))/energy_normalization)
%disp(min(Mu_pathGM(:, 1)/energy_normalization))
%disp(max(Mu_pathGM(:, 1)/energy_normalization))

% figure(4)
% mesh(xlat_rsforU, ylat_rsforU, abs(Wannier_2D_Norm))
% hold on
% mesh(xlat_rsforU, ylat_rsforU, abs(Wannier_a1))
% hold on
% mesh(xlat_rsforU, ylat_rsforU, abs(Wannier_a2))
% hold on
% mesh(xlat_rsforU, ylat_rsforU, abs(Wannier_a3))
% hold off
% xlim([-200,200])
% ylim([-200,200])
% axis('square')
% 
% disp(['t from BZ = {' num2str(t2_a1) ', '  num2str(t2_a2) '}'])
% disp(['t from space int by FD = {' num2str(t_a1) ', '  num2str(t_a2) '}'])
% disp(['t from space int by DFT = {' num2str(t_ALT_a1) ', '  num2str(t_ALT_a2) '}'])
% disp(['t from GMK spaghetti = ' num2str(t_fromGMKspaghetti)])
% disp(['V for tight binding approximation = ' num2str(T_a1)])

% outtt = [b_array' T_a1_energy_array' V_graphene_energy_array' V_He_He_energy_array' V_old_array' V_sym_array'];
% 
% writematrix(real(outtt), 'b_energy_V.txt')

% xxx = ([2.574 :0.001: 3.965] - 3.965); 
% figure(114); 
% plot(L_pathKG, ...
%      (Mu_pathKG(:, 1)+0.8867)*5/3.148 );
% hold on
% plot(1.72*xxx+ 1.76 * 3.965, (-(cos(xxx/3.3208e-01) + 2*cos(xxx/2/3.3208e-01)) + 3)*...
%     1.1/4.5 - 110.9, 'g'  )
 
% 
% figure(114);
% mesh(x/Lx_singlePeriod, y/Ly_singlePeriod, Wannier_combined);
% xlabel('x'); ylabel('y'); zlabel('1D Wannier combined')


% Plot the difference between the moduli of the direct 2D and "synthesized" 2D Wanniers:
% figure(115);
% mesh(x/Lx_singlePeriod, y/Ly_singlePeriod, Wannier_combined - abs(Wannier_2D_Norm));
% xlabel('x'); ylabel('y'); zlabel('|"1D" Wannier| - |Wannier|')

% ----------------------------------------------------------------------------------------------
% ----------------------------------------------------------------------------------------------
% Changelog:
% v4: 
%    Added more periods in x and y for calculation of Wannier function
% v5: 
%    1) Corrected a minor error in OuterExp. 
%    2) Wannier function is obtained by integration over whole Brillouin zone.
%    3) Wavenumbers (denoted here by k1, k2) in the Brillouin zone are computed
%       *not* in [-edge, +edge-dk], as would be usual for Fourier transform,
%       but in [-edge+dk/2, +edge-dk/2].
%       This is done to avoid the situation where at the "-edge" the Hamiltonian matrix
%       would have double eigenvalues and hence there would be an ambiguity of which of
%       the two corresponding eigenvectors one would need to use for the Bloch function.
%       (This is explained in more detail in the appropriate place in the code.)
%    4) The earlier problem related to the fact that Matlab computes an eigenvector only
%       up to the sign, is now fixed. (This is described in detail in the appropriate 
%       place in the code.) Before this issue was fixed, the Wannier function would look
%       far from correct. 
%       NOTE / OPEN ISSUE (not impacting this code):
%       The Hamiltonian matrix is real symmetric, and hence its eigenvectors are real.
%       Then, since Matlab normalizes them, the only ambiguity is in their sign, which
%       was easy to fix. However, if we ever have to deal with a complex Hermitian
%       Hamiltonian matrix, the corresponding eigenvectors can be complex. Then the 
%       ambiguity can be in the phase (which is now continuous, rather than just 0 or pi)
%       of the vector. This is the aforementioned open issue. 
% v6:
%    1) Added an alternative way of creating the Hamiltonian matrix. 
%       This way is more intuitive to set up, as it follows *directly* from the equation
%       for the Fourier coefficients of the Bloch function.
%       However, it does require us to adjust the way Matlab reshapes matrices with the
%       intuitive way in which we set up our matrix. Fortunately, this can be done with
%       just one Matlab command. Details inside. 
%    2) Added the capability where we can compare our result with the product of 1D results
%       in Adrian's code. 
% v6.5:
%    1) Fixed the way we defined [Kn, Km] = meshgrid(Knvec, Kmvec). This is the matrix of
%       dimension (2Nx+1) x (2Ny+1) (see below). 
%       In v6 and earlier, it was defined as [Km, Kn] = meshgrid(Kmvec, Knvec), which has
%       a wrong dimension. As a result, in v6 we had to permute dimensions of HM before
%       reshaping it; otherwise, results would be off. With the above fix of [Kn,Km],
%       the plain  reshape  with only one additional trick - the transposition (see the
%       EXPLANATION in the code) correctly flattens  HM. 
%    2) Introduced the ability to have different numbers of Fourier harmonics inx- and y-
%       directions. This was not so much from the desire to have those numbers to be different,
%       but rather to make sure that we correctly treat x- and y-dimensions. 
%  v7:
%     1) Now loads the graphene potential from txt file
%  v8:
%     1) Changes the size of the potential grid loaded from file, now the
%        same size as one primitive cell.
%  v9:
%     1) Change of variables from x-y to along the lattice vectors a1 and
%     a2
%     2) Loads fourier terms, potential, and lattice directly from file
%  v9a:
%     Test case for cos potential
%  v10:
%     1) Adds calculation of U, the on-site interaction
%  v11:
%     Includes calculation of V, nearest neighbor interactions.
%  v13:
%     Only Calculation of t
