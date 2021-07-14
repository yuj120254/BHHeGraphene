% graphene_lattice class
% contains all the parameters and lattice vectors of strained graphene
classdef Graphene_Lattice
    properties
        a0 % distance between carbon atoms in graphene
        sigma  % potential fitting parameter
        delta  % Strain
        a1 % Lattice Vectors
        a2 %
        g1 % Reciprocal Lattice Vectors
        g2 %
        A % area of Brillouin zone
        b1 % Basis vectors
        b2 %
    end
    
    methods 
        function glat = Graphene_Lattice(sigma_temp, delta_temp, a0_temp, a1_temp, a2_temp)
            if nargin == 5
                glat.a0 = a0_temp;
                glat.sigma = sigma_temp;
                glat.delta = delta_temp;
                glat.a1 = a1_temp;
                glat.a2 = a2_temp;
                glat.g1 = [a2_temp(2)/(a1_temp(1)*a2_temp(2)-a1_temp(2)*a2_temp(1)), a2_temp(1)/(a1_temp(2)*a2_temp(1)-a1_temp(1)*a2_temp(2))].*(2*pi);
                glat.g2 = [a1_temp(2)/(a1_temp(2)*a2_temp(1)-a1_temp(1)*a2_temp(2)), a1_temp(1)/(a1_temp(1)*a2_temp(2)-a1_temp(2)*a2_temp(1))].*(2*pi);
                glat.A = sqrt(sum(cross([a1_temp, 0], [a2_temp, 0]).^2));
                glat.b1 = [a2_temp(1), -1/2*(1 + delta_temp)*a0_temp];
                glat.b2 = [0, -a0_temp*(1 + delta_temp)];
            end
        end
        
        function F = f(k1, k2, glat)
            F = exp(1i*(dot(glat.g1*k1 + glat.g2*k2, glat.b1))) + exp(1i*(dot(glat.g1*k1 + glat.g2*k2, glat.b2)));
        end

        function G = g(k1, k2, glat)
            G = sqrt(sum((k1*glat.g1 + k2*glat.g2).^2));
        end

        function c = C(k1, k2, z, glat)
            gC = g(k1, k2, glat);

            lj_factor = gC*glat.sigma^2/(2*z);

            c = 1/30*(lj_factor)^5*besselk(5, gC*z) - 2*(lj_factor)^2*besselk(2, gC*z);
        end

        function locexp = location_exponent(k1, k2, r, glat)
            locexp = exp(1i*(dot(glat.g1.*k1 + glat.g2.*k2, r)));
        end

        function un = Un(k1, k2, x, y, z, glat)
            r = [x, y];
            e_gr = location_exponent(k1, k2, r, glat);

            un =  f(k1, k2, glat) * C(k1, k2, z, glat) * e_gr;
        end

        function vn = Vn(x, y, z, gterms, glat)
            V = 0;

            for i = -gterms:gterms
                for j = -gterms:gterms
                    if (i == 0) && (j == 0)
                        V = V + C0(z, glat);
                    else
                        V = V + Un(i, j, x, y, z, glat);
                    end
                end
            end

            vn = V;
        end	

        function c0 = C0(z, glat)
            c0 = 4/5*(glat.sigma/z)^10 - 2*(glat.sigma/z)^4;
        end

        function vz = V_z_0(z_nondim, gterms, glat)
            V = zeros(length(z_nondim), 1);

            for i = -gterms:gterms
                for j = -gterms:gterms
                    for iz = 1:length(z_nondim)
                        if (i == 0) && (j == 0)
                            V(iz) = V(iz) + C0(z_nondim(iz), glat);
                        else
                            V(iz) = V(iz) + f(i,j, glat)*C(i,j,z_nondim(iz), glat);
                        end
                    end
                end
            end

            vz = V;
        end

        function vr = V_r(x, y, z, gterms, glat)
            V = zeros(length(x), length(y));

            for i = 1:length(x)
                for j = 1:length(y)
                    V(i,j) = real(Vn(x(i), y(j), z, gterms, glat));
                end
            end

            vr = V;
        end

        function F_t = Fourier_Terms(gterms, zmin, glat)
            Vfourier = complex(zeros(2*gterms+1, 2*gterms+1));

            for k1 = -gterms:gterms
                for k2 = -gterms:gterms
                    if (k1 == 0 && k2 == 0)
                        Vfourier(k1+gterms+1,k2+gterms+1) = C0(zmin, glat);
                    else
                        Vfourier(k1+gterms+1,k2+gterms+1) = f(k1, k2, glat) * C(k1, k2, zmin, glat);
                    end
                end
            end

            F_t = Vfourier;
        end

        function v = V_fourier(x, y, gterms, Vfourier, scale, glat)
            V = complex(zeros(length(x), length(y)));

            for i = 1:length(x)
                for j = 1:length(y)
                    r = [x(i), y(j)];
                    for k1 = -gterms:gterms
                        for k2 = -gterms:gterms
                            e_gr = location_exponent(k1, k2, r, glat);
                            V(i,j) = V(i,j) + scale*Vfourier(k1+gterms+1,k2+gterms+1)*e_gr;
                        end
                    end
                end
            end

            v = V;
        end	

        function v = V_fourier_grid(xgrid, ygrid, gterms, Vfourier, scale, glat)
            V = complex(zeros(size(xgrid, 1), size(ygrid, 2)));

            for ix = 1:size(xgrid, 1)
                for iy = 1:size(ygrid, 2)
                    r = [xgrid(ix, iy), ygrid(ix, iy)];
                    for k1 = -gterms:gterms
                        for k2 = -gterms:gterms
                            location_exponent = exp(1i*(dot(glat.g1.*k1 + glat.g2.*k2, r)));
                            V(ix,iy) = V(ix,iy) + scale*Vfourier(k1+gterms+1,k2+gterms+1)*location_exponent;
                        end
                    end
                end
            end

            v = V;
        end	
    end
end

