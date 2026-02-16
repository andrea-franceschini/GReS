function [E_rock, nu_rock, azimuths, dips, Kn, Ks, numbers_of_joints] = ...
    set_parameters(rock, jointFamilies)

    % Rock properties
    E_rock = rock.E;
    nu_rock = rock.nu;

    keys = fieldnames(jointFamilies);
    nFam = length(keys);

    azimuths.mean = zeros(1,nFam);
    azimuths.std  = zeros(1,nFam);
    dips.mean     = zeros(1,nFam);
    dips.std      = zeros(1,nFam);

    Kn.mean = zeros(1,nFam);
    Kn.std  = zeros(1,nFam);
    Kn.min  = zeros(1,nFam);
    Kn.max  = zeros(1,nFam);

    Ks.mean = zeros(1,nFam);
    Ks.std  = zeros(1,nFam);
    Ks.min  = zeros(1,nFam);
    Ks.max  = zeros(1,nFam);

    for i = 1:nFam
        jf = jointFamilies.(keys{i});

        azimuths.mean(i) = jf.azimuth_mean;
        azimuths.std(i)  = jf.azimuth_std;

        dips.mean(i) = jf.dip_mean;
        dips.std(i)  = jf.dip_std;

        Kn.mean(i) = jf.Kn_mean;
        Kn.std(i)  = jf.Kn_std;
        Kn.min(i)  = jf.Kn_min;
        Kn.max(i)  = jf.Kn_max;

        Ks.mean(i) = jf.Ks_mean;
        Ks.std(i)  = jf.Ks_std;
        Ks.min(i)  = jf.Ks_min;
        Ks.max(i)  = jf.Ks_max;
    end

    % Number of joints (CSV string -> numeric array)
    first = jointFamilies.(keys{1}).number_of_joints;
    numbers_of_joints = zeros(length(first),nFam);
    numbers_of_joints(:,1) = first;
    for j = 2 : nFam
        col = jointFamilies.(keys{j}).number_of_joints;
        numbers_of_joints(:,j) = col;
    end

end
