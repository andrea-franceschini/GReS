cd(fullfile('clipper2','private'));
% original - mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp')
mex('clipper.cpp','mexclipper.cpp');