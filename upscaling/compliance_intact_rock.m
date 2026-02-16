function S = compliance_intact_rock(E, nu)

    G = E / (2*(1+nu));
    S = zeros(6);

    S(1,1) = 1/E;
    S(2,2) = 1/E;
    S(3,3) = 1/E;

    S(1,2) = -nu/E; S(1,3) = -nu/E;
    S(2,1) = -nu/E; S(2,3) = -nu/E;
    S(3,1) = -nu/E; S(3,2) = -nu/E;

    S(4,4) = 1/G;
    S(5,5) = 1/G;
    S(6,6) = 1/G;
end
