% compile VTKread library using CMAKE
% if ~exist('build','dir')
   mkdir build
   cd build
   system('cmake ..');
   system('make');
%end