function [mesh, rock, blockSize, jointFamilies, druckerPrager, ...
    inputFiles] = readInputFile(fileName)

    problem = readstruct(fileName,AttributeSuffix="");

    mesh = problem.Mesh;
    rock = problem.Rock;
    blockSize = problem.JointFamilies.BlockSize;
    jointFamilies = problem.JointFamilies.Families;
    druckerPrager = problem.DruckerPrager;
    inputFiles = problem.InputFiles;

    fields = fieldnames(jointFamilies);
    for idx = 1 : length(fields)
        jointFamilies.(fields{idx}).number_of_joints = ...
            str2num(jointFamilies.(fields{idx}).number_of_joints);
    end

end
