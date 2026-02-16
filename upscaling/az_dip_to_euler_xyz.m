function [alpha, beta, gamma] = az_dip_to_euler_xyz(azimuth_deg, dip_deg)

    A = deg2rad(90.0 - azimuth_deg);
    D = deg2rad(dip_deg);

    Rx = [1 0 0;
          0 cos(D) -sin(D);
          0 sin(D)  cos(D)];

    Rz = [cos(A) -sin(A) 0;
          sin(A)  cos(A) 0;
          0 0 1];

    R = Rz * Rx;

    beta  = asin(-R(1,3));
    alpha = atan2(R(3,2), R(3,3));
    gamma = atan2(R(2,1), R(1,1));

    alpha = rad2deg(alpha);
    beta  = rad2deg(beta);
    gamma = rad2deg(gamma);
end
