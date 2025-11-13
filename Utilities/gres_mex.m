function gres_mex(mxfName, varargin)
% compileMex - Compiles a MEX file if it does not already exist.
%
% Usage:
%   compileMex(fName)
%   compileMex(fName, cppFile1, cppFile2, ...)
%
% Inputs:
%   mxfName     - (string or char) base name of the MEX interface file
%   varargin  - additional .cpp files to include in compilation

    % Ensure fName is a string
    mxfName = string(mxfName);

    % Check if the MEX file already exists
    if exist(mxfName, "file") == 3
        fprintf("Mex file '%s' already exists. Compilation is skipped.\n", mxfName);
        return;
    end

    % Collect all source files
    cppFiles = [{mxfName + ".cpp"}, varargin];

    mexArgs = ["-O", cppFiles{:}];
    mex(mexArgs{:});

end