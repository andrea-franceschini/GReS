function grdecl = movePillars(grdecl, dir, target, geom_fac, varargin)

% input is a corner point grid as obtained from mrst

% move the pillars according to a geometric distribution refinement around
% SOME target pillars

% for each set of pillars along direction x, split the set around the
% target 

% determine the modified pillars length compare to the initial ones.
% translate the pillars accordingly

% use geometric factor to determine the distribution

% varargin: string array specifying where to refine wrt the target (before, after, both)
if isempty(varargin)
  side = ["before", "after"];
else
  side = varargin{1};
end

nPillarsDir = grdecl.cartDims(dir)+1;
if dir == 1 
  dirNorm = 2;
else
  dirNorm = 1;
end

nPillarsNorm = grdecl.cartDims(dirNorm)+1;

% reshape coord depending on the refinement direction
if dir == 2
  N = nPillarsNorm;
  c = reshape(grdecl.COORD,6,N,[]);
  c = permute(c,[1,3,2]);
  %c = c(:);
elseif dir == 1
  N = nPillarsDir;
  c = reshape(grdecl.COORD,6,N,[]);
end




for in = 1:nPillarsNorm

  % get total lenght
  c1 = c(:,1,in);
  c2 = c(:,target,in);
  cEnd = c(:,end,in);

  %average distance between first and target pillar
  L_left = mean(abs(c1(dir:3:end)-c2(dir:3:end)));
  L_right = mean(abs(cEnd(dir:3:end)-c2(dir:3:end)));

  Nseg_left = target-1;
  Nseg_right = nPillarsDir - target;
  L1_left = (L_left * (1-geom_fac))/(1 - geom_fac^(Nseg_left));
  L1_right = (L_right * (1-geom_fac))/(1 - geom_fac^(Nseg_right));
  lengths_left = L1_left*geom_fac.^(0:Nseg_left-1);
  lengths_right = flip(L1_right*geom_fac.^(0:Nseg_right-1));

  if ismember("before",side)
    for i=2:target-1
      % move pillars to match target lenghts
      d = mean(abs( c(dir:3:end,i,in)-c1(dir:3:end)));
      c(dir:3:end,i,in) = c(dir:3:end,i,in) - d + sum(lengths_left(1:i-1));
    end
  end

  if ismember("after", side)
    for i=target+1:nPillarsDir-1
      d = mean(abs( c(dir:3:end,i,in)-c2(dir:3:end)));
      c(dir:3:end,i,in) = c(dir:3:end,i,in) - d + sum(lengths_right(1:i-target));
    end
  end


end

% permute back
if dir == 2
  c = permute(c,[1,3,2]);
end

% redefine coord

grdecl.COORD = c(:);

end








