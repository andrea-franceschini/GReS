function log = gres_log()
log = getappdata(0,'gres_log');
assert(~isempty(log),['GReS logger not defined. Run initGReS to initialize ' ...
  'the simulator'])
end
