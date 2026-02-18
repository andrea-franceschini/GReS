function [E_rock, nu_rock, azimuths, dips, Kn, Ks, numbers_of_joints] = ...
    set_parameters(rock, jointFamilies)

    % Rock properties
    E_rock = rock.E;
    nu_rock = rock.nu;

    keys = fieldnames(jointFamilies);
    nFam = length(keys);

    azimuths.mean = zeros(nFam,1);
    azimuths.std  = zeros(nFam,1);
    dips.mean     = zeros(nFam,1);
    dips.std      = zeros(nFam,1);

    Kn.mean = zeros(nFam,1);
    Kn.std  = zeros(nFam,1);
    Kn.min  = zeros(nFam,1);
    Kn.max  = zeros(nFam,1);

    Ks.mean = zeros(nFam,1);
    Ks.std  = zeros(nFam,1);
    Ks.min  = zeros(nFam,1);
    Ks.max  = zeros(nFam,1);

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
