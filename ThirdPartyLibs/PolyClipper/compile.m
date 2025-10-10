cd(fullfile('clipper2','private'));
mex('-D__int64=__int64_t','clipper.cpp','mexclipper.cpp')
