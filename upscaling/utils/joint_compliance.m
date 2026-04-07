function S = joint_compliance(Kn, Ks)

    S = [1/Ks 0 0;
         0 1/Ks 0;
         0 0 1/Kn];
end
