function root = gres_root()
root = getappdata(0,'gres_root');
assert(~isempty(root),['GReS root not defined. Run initGReS to initialize ' ...
  'the simulator'])
end
