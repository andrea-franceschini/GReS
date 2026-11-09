bounds = [0 1; 0 1; 0 1];

% Persistent steep fracture at x = 0.48.
main = [.48 .12 .08;
        .48 .88 .08;
        .48 .88 .92;
        .48 .12 .92];

% Cross-cutting steep fracture at y = 0.57.
crossCut = [.10 .57 .20;
            .86 .57 .20;
            .86 .57 .82;
            .10 .57 .82];

% Interior branch that terminates against the main fracture.
branch = [.17 .32 .28;
          .48 .32 .28;
          .48 .32 .72;
          .17 .32 .72];

% Inclined connector. It cuts the main fracture and the interior branch.
inclined = [.28 .24 .25;
            .76 .24 .25;
            .76 .51 .79;
            .28 .51 .79];

% Short parallel fracture, offset 0.06 from the main plane.
parallel = [.54 .19 .18;
            .54 .44 .18;
            .54 .44 .62;
            .54 .19 .62];

% Branch that reaches the box boundary and terminates at the main plane.
boundaryBranch = [0   .76 .36;
                  .48 .76 .36;
                  .48 .76 .67;
                  0   .76 .67];

network = {main, crossCut, branch, inclined, parallel, boundaryBranch};
names = ["Matrix", "Main", "Cross-cut", "Interior branch", ...
         "Inclined connector", "Parallel", "Boundary branch"];

[matrix, fractures] = createFracturedMesh(bounds, network, ...
    'matrixMeshSize', 0.10, 'fractureMeshSize', 0.025, ...
    'transitionFactor', 1.0);

plotMesh({matrix, fractures{:}}, 'damage_zone', names);

   