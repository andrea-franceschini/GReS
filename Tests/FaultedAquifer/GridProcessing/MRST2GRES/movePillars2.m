function grdecl = movePillars2(grdecl, dir, targets, widths, hmin, hmax)

% input is a corner point grid as obtained from mrst

% move the pillars according to a geometric distribution refinement around
% SOME target pillars

% for each set of pillars along direction x, split the set around the
% target 

% determine the modified pillars length compare to the initial ones.
% translate the pillars accordingly

% use geometric factor to determine the distribution

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
  cEnd = c(:,end,in);

  %average distance between first and target pillar
  L = mean(abs(c1(dir:3:end)-cEnd(dir:3:end)));

  lengths_seg = multiBumpSpacings(nPillarsDir, L, targets, widths, hmin, hmax);

  for i=2:nPillarsDir-1
    % move pillars to match target lenghts
    d = mean(abs( c(dir:3:end,i,in)-c1(dir:3:end)));
    c(dir:3:end,i,in) = c(dir:3:end,i,in) - d + sum(lengths_seg(1:i-1));
  end

end


% permute back
if dir == 2
  c = permute(c,[1,3,2]);
end

% redefine coord

grdecl.COORD = c(:);

end



function dx = multiBumpSpacings(N, L, targets, widths, hmin, hmax)
% Returns N-1 segment lengths along a line of length L with multiple refinement bumps
% targets, widths: normalized positions and widths in [0,1]

% Sample finely
s = linspace(0,1,1000);
h = hmax*ones(size(s));
for i = 1:length(targets)
%     h = min(h, hmax - (hmax-hmin)*exp(-((s-targets(i))/widths(i)).^2/2));
%      h = min(h, hmin + (hmax-hmin)*exp(-((s-targets(i))/widths(i)).^2/2));
 h = min(h, hmax - (hmax-hmin)*exp(-((s-targets(i))/widths(i)).^2/2));
end

% Compute cumulative distribution
cs = cumsum(h);
cs = cs / cs(end);  % normalize to [0,1]

% Interpolate node positions at uniform fractions
nodes = interp1(cs, s, linspace(0,1,N), 'linear');

nodes(1) = 0;

% Compute spacings
dx = diff(nodes) * L;
end




