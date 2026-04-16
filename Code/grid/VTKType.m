classdef VTKType < uint8

  % VTKTypes implemented in GReS

  enumeration
    Tri (5)           % 3 node triangle
    Polygon (7)       % n-node polygon
    Quad  (9)         % 4 node quadrilateral
    Quad9 (28)        % 9 node quadrilateral (biquadratic)
    Tetra (10)        % 4 node tetrahedron
    Hexa  (12)        % 8 node hexahedron
    Hexa27 (29)       % 27 node hexahedron (triquadratic)


    % LIST OF VTK TYPES FOR REFERENE
    %  1  VTK_VERTEX
    %  2  VTK_POLY_VERTEX
    %  3  VTK_LINE
    %  4  VTK_POLY_LINE
    %  5  VTK_TRIANGLE
    %  6  VTK_TRIANGLE_STRIP
    %  7  VTK_POLYGON
    %  8  VTK_PIXEL
    %  9  VTK_QUAD
    % 10  VTK_TETRA
    % 11  VTK_VOXEL
    % 12  VTK_HEXAHEDRON
    % 13  VTK_WEDGE
    % 14  VTK_PYRAMID
    % 28  VTK_BIQUADRATIC_QUAD
    % 29  VTK_TRIQUADRATIC_HEXAHEDRON
  end

  methods (Static)


    function t2d = to2D(t3d)
      % Convert 3D VTK type(s) to corresponding face type(s)

      v = uint8(t3d);

      t2d = zeros(size(v), 'uint8');

      mask = (v == uint8(VTKType.Tetra));
      t2d(mask) = uint8(VTKType.Tri);

      mask = (v == uint8(VTKType.Hexa));
      t2d(mask) = uint8(VTKType.Quad);

      mask = (v == uint8(VTKType.Hexa27));
      t2d(mask) = uint8(VTKType.Quad9);

      % Check for invalid inputs
      if any(t2d == 0)
        invalid = unique(v(t2d == 0));
        error("Unsupported 3D VTK type(s): %s", mat2str(invalid));
      end

      % Convert back to enum
      t2d = VTKType(t2d);
    end



    function t3d = to3D(t2d)
      % Convert 3D VTK type(s) to corresponding face type(s)

      v = uint8(t2d);

      t3d = zeros(size(v), 'uint8');

      mask = (v == uint8(VTKType.Tri));
      t3d(mask) = uint8(VTKType.Tetra);

      mask = (v == uint8(VTKType.Quad));
      t3d(mask) = uint8(VTKType.Hexa);

      mask = (v == uint8(VTKType.Quad9));
      t3d(mask) = uint8(VTKType.Hexa27);

      % Check for invalid inputs
      if any(t3d == 0)
        invalid = unique(v(t3d == 0));
        error("Unsupported 2D VTK type(s): %s", mat2str(invalid));
      end

      % Convert back to enum
      t3d = VTKType(t3d);
    end



    function t3d = types3D()
      t3d = [10,12,29];
    end

    function t2d = types2D()
      t2d = [5,9,28];
    end

  end

end
