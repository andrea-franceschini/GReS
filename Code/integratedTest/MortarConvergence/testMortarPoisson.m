% study parameters
for elem_type = ["hexa","hexa27"]

  for integration_type = ["SegmentBasedQuadrature",...
                          "RBFquadrature",...
                          "ElementBasedQuadrature",...
                          ]

    runConvPoisson;

    msg = "Error for %s element with %s mortar quadrature scheme";

    switch elem_type
      % ensure correct convergence rate 
      case "hexa"
        assert(all([L2ord>1.8;L2ord<2.5]),...
          msg,elem_type,integration_type)
        assert(all([H1ord>0.8;H1ord<1.5]),...
          msg,elem_type,integration_type)
      case "hexa27"
        assert(all([L2ord>2.8;L2ord<3.5]),...
          msg,elem_type,integration_type)
        assert(all([H1ord>1.8;H1ord<2.5]),...
          msg,elem_type,integration_type)
    end
  end
  clearvars
end
