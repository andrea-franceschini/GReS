function log = gresLog()
log = getappdata(0,'gresLog');
assert(~isempty(log),['GReS logger not defined. Run initGReS to initialize ' ...
  'the simulator'])
end
