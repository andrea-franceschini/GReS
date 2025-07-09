clear
close all
clc
%
h = 1./[3 6 12 15 18 21 24];
integration_type = ["SegmentBased", "RBF", "ElementBased"];
figure(1)

lw = 1;
ms = 8;


for i_t = integration_type
  switch i_t
    case 'SegmentBased'
      col = 'k';
      m = 'o';
      ls = '-';
      lw = 2;
      nInt = 0;
    case 'ElementBased'
      col = 'r';
      m = 's';
      ls = '--';
      lw = 1;
      nInt = 0;
    case 'RBF'
      col = 'b';
      m = '+';
      ls = '-';
      lw = 1;
      nInt = [4 5 6];
  end

  if strcmp(i_t,'SegmentBased')
    nG = 7;
  else
    nG = [4 6 16];
  end
  for ngp = nG
    for n_i = nInt
      fname = i_t+"_"+num2str(ngp);
      if i_t=="RBF"
        fname = fname + "_" + num2str(n_i);
      end
      fname = fname+".mat";
      out = load(fullfile("Output/OUT_HEXA27_2",fname));
      L2 = out.L2norm;
      H1 = out.H1norm;
      nref = numel(L2);
      leg_txt = i_t+" - "+num2str(ngp)+" GP";
      if i_t == "RBF"
        leg_txt = leg_txt + " - " + num2str(n_i) + " interp points";
      end
      loglog(h(1:nref),L2,...
        "DisplayName",leg_txt,...
        'LineStyle',ls,...
        'Color',col,...
        'Marker',m,...
        'LineWidth',lw,...
        'MarkerSize',ms,...
        'MarkerEdgeColor',col)
      hold on
      loglog(h(1:nref),H1,...
        "DisplayName",leg_txt,...
        'LineStyle',ls,...
        'Color',col,...
        'Marker',m,...
        'LineWidth',lw,...
        'MarkerSize',ms,...
        'MarkerEdgeColor',col)
    end
  end
end