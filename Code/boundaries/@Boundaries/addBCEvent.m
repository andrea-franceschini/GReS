function addBCEvent(boundaries,name,varargin)
% ADDBCEVENT  Add a time-dependent event to an existing boundary condition.
%
%   ADDBCEVENT(obj, name, Key, Value) appends a (time, value) pair to the
%   boundary condition object boundaries. Events are kept sorted in
%   ascending time order. Duplicate time entries are not allowed.
%
% -------------------------------------------------------------------------
% INPUT ARGUMENTS (Name-Value Pairs)
% -------------------------------------------------------------------------
%
%   'time'   - (numeric) Time at which the BC value is applied.
%              If time is a function handle f(t), it represent a time
%              varying scalar factor multiplying the associated value.

%
%   'value'  - BC value associated with the given time. 
%              Can be a single scalar, a list of scalars mathcing the BC
%              entities, a space dependent function f(x,y,z) or a path to a
%              file containing the list of values. 
%
% -------------------------------------------------------------------------
% NOTES
% -------------------------------------------------------------------------
%
%   - Events are automatically sorted by ascending time after each call. -
%   - Registering two events with the same time raises an error. 
%
% -------------------------------------------------------------------------
% SEE ALSO
% -------------------------------------------------------------------------
%   addBC

boundaries.getData(name).data.addBCEvent(varargin{:});

end

