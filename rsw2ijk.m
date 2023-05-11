function [r, v, u] = rsw2ijk(r_rsw, v_rsw, omega, i, w, f)

    function mat1 = rot1(x)
        mat1 = [1    0      0;
               0  cos(x) sin(x);
               0 -sin(x) cos(x)];
    end

    function mat3 = rot3(x)
        mat3 = [cos(x) sin(x) 0;
               -sin(x) cos(x) 0;
                   0       0  1];
    end

u = f + w;
rot313 = rot3(-omega) * rot1(-i) * rot3(-u);

r = rot313 * r_rsw;
v = rot313 * v_rsw;

end



