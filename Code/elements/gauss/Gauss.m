classdef Gauss < handle
  %UNTITLED Summary of this class goes here
  %   Detailed explanation goes here

  properties (Access = public)
    coord;
    weight;
    nNode;
  end

  methods (Access = public)
    function obj = Gauss(cType,nG)
      %UNTITLED Construct an instance of this class
      %   Detailed explanation goes here
      if nargin ~= 2
        error('Incorrect enough input argument in Gauss class');
      end
      setGaussPoints(obj,cType,nG);
    end
  end

  methods (Access = private)

    function setGaussPoints(obj,cType,nG)
      %METHOD1 Summary of this method goes here
      %   Detailed explanation goes here
      switch cType
        case 5 % Triangle
          [c,w] = Gauss.pointsTriangle(nG);
          obj.coord = c;
          obj.weight = w;
          obj.nNode = numel(w);
        case {9,28} % Quadrilateral
          [c1,w1] = Gauss.points1D(nG);
          [y, x] = meshgrid(c1,c1);
          obj.coord = [x(:), y(:)];
          obj.weight = bsxfun(@(a,b) a*b,w1,(w1)');
          obj.weight = obj.weight(:);
          obj.nNode = (numel(w1))^2;
        case 10 % Tetrahedron
          [c,w] = Gauss.pointsTetra(nG);
          obj.coord = c;
          obj.weight = w;
          obj.nNode = numel(w);
        case {12,29} % Hexahedron
          [c1,w1] = Gauss.points1D(nG);
          [y, x, z] = meshgrid(c1,c1,c1);
          obj.coord = [x(:), y(:), z(:)];
          weightGaussTmp = bsxfun(@(a,b) a*b,w1,(w1)');
          obj.weight = bsxfun(@(a,b) a*b,weightGaussTmp(:),(w1)');
          obj.weight = obj.weight(:);
          obj.nNode = (numel(w1))^3;
        otherwise
          error('Gauss points for VTK element type %d not yet implemented /n',cType);
      end
    end
  end


  methods (Static)

   function [coord1D, weight1D] = points1D(nG)

  % get 1D points location for tensor product rule
  switch nG
    case 1
      coord1D = 0.0;
      weight1D = 2.0;
    case 2
      l = 1/sqrt(3);
      coord1D = [-l; l];
      weight1D = ones(2,1);
    case 3
      a = sqrt(0.6);
      b = 0.0;
      wa = 5/9;
      wb = 8/9;
      coord1D = [-a; b; a];
      weight1D = [wa; wb; wa];
    case 4
      a = 0.861136311594953;
      b = 0.339981043584856;
      wa = 0.347854845137454;
      wb = 0.652145154862546;
      coord1D = [-a; -b; b; a];
      weight1D = [wa; wb; wb; wa];
    case 5
      a = 0.906179845938664;
      b = 0.5384693101056831;
      c = 0.0;
      wa = 0.2369268850561891;
      wb = 0.4786286704993665;
      wc = 0.568888888888888;
      coord1D = [-a; -b; c; b; a];
      weight1D = [wa; wb; wc; wb; wa];
    case 6
      a = 0.9324695142031521;
      b = 0.6612093864662645;
      c = 0.2386191860831969;
      wa = 0.1713244923791704;
      wb = 0.3607615730481386;
      wc = 0.4679139345726910;
      coord1D = [-a; -b; -c; c; b; a];
      weight1D = [wa; wb; wc; wc; wb; wa];
    case 8
      a = 0.1834346424956498;
      b = 0.5255324099163290;
      c = 0.7966664774136267;
      d = 0.9602898564975363;
      wa = 0.3626837833783620;
      wb = 0.3137066458778873;
      wc = 0.2223810344533745;
      wd = 0.1012285362903763;
      coord1D = [-d; -c; -b; -a; a; b; c; d];
      weight1D = [wd; wc; wb; wa; wa; wb; wc; wd];
    case 16
      a = 0.9894009349916499;
      b = 0.9445750230732326;
      c = 0.8656312023878318;
      d = 0.7554044083550030;
      e = 0.6178762444026438;
      f = 0.4580167776572274;
      g = 0.2816035507792589;
      h = 0.0950125098376374;
      wa = 0.0271524594117541;
      wb = 0.0622535239386479;
      wc = 0.0951585116824928;
      wd = 0.1246289712555339;
      we = 0.1495959888165767;
      wf = 0.1691565193950025;
      wg = 0.1826034150449236;
      wh = 0.1894506104550685;
      coord1D = [-a; -b; -c; -d; -e; -f; -g; -h; h; g; f; e; d; c; b; a];
      weight1D = [wa; wb; wc; wd; we; wf; wg; wh; wh; wg; wf; we; wd; wc; wb; wa];
    otherwise
      error(['Gauss rule for %i numb. of GP is not available. Available ' ...
        'rules for 1,2,3,4,5,6,8,16 GP\n'],nG);
  end
end
    function [coord,weight] = pointsTriangle(nG)

      % get 1D points location for reference triangle
      % from Dunavant, 1984
      parentArea = 0.5;
      switch nG
        case 1  % Degree 1 exact, 1 point
          coord = [1/3, 1/3];
          weight = 1.0;

        case 3  % Degree 2 exact, 3 points
          coord = [
            0.166666666666667, 0.166666666666667;
            0.666666666666667, 0.166666666666667;
            0.166666666666667, 0.666666666666667
            ];
          weight = repmat(1/3, 3, 1);

        case 4  % Degree 3 exact, 4 points
          coord = [
            1/3, 1/3;
            0.6, 0.2;
            0.2, 0.6;
            0.2, 0.2
            ];
          weight = [
            -27/48;
            25/48;
            25/48;
            25/48
            ];

        case 6  % Degree 4 exact, 6 points
          coord = [
            0.445948490915965, 0.445948490915965;
            0.445948490915965, 0.108103018168070;
            0.108103018168070, 0.445948490915965;
            0.091576213509771, 0.091576213509771;
            0.091576213509771, 0.816847572980459;
            0.816847572980459, 0.091576213509771
            ];
          weight = [
            0.223381589678011;
            0.223381589678011;
            0.223381589678011;
            0.109951743655322;
            0.109951743655322;
            0.109951743655322
            ];

        case 7  % Degree 5 exact, 7 points
          coord = [
            1/3, 1/3;
            0.470142064105115, 0.470142064105115;
            0.470142064105115, 0.059715871789770;
            0.059715871789770, 0.470142064105115;
            0.101286507323456, 0.101286507323456;
            0.101286507323456, 0.797426985353087;
            0.797426985353087, 0.101286507323456
            ];
          weight = [
            0.225000000000000;
            0.132394152788506;
            0.132394152788506;
            0.132394152788506;
            0.125939180544827;
            0.125939180544827;
            0.125939180544827
            ];

        case 12  % Degree 6 exact, 12 points
          coord = [
            0.249286745170910, 0.249286745170910;
            0.249286745170910, 0.501426509658179;
            0.501426509658179, 0.249286745170910;

            0.063089014491502, 0.063089014491502;
            0.063089014491502, 0.873821971016996;
            0.873821971016996, 0.063089014491502;

            0.310352451033785, 0.636502499121399;
            0.636502499121399, 0.053145049844816;
            0.053145049844816, 0.310352451033785;
            0.636502499121399, 0.310352451033785;
            0.310352451033785, 0.053145049844816;
            0.053145049844816, 0.636502499121399
            ];

          weight = [
            0.116786275726379;
            0.116786275726379;
            0.116786275726379;

            0.050844906370207;
            0.050844906370207;
            0.050844906370207;

            0.082851075618374;
            0.082851075618374;
            0.082851075618374;
            0.082851075618374;
            0.082851075618374;
            0.082851075618374
            ];
        otherwise
          error(['Unsupported number of gauss points for triangle quadrature.' ...
            'Available integration order are: \n' ...
            ' order 1: 1 gp \n order 2: 3 gp \n order 3: 4 gp \n' ...
            ' order 4: 6 gp \n order 5: 7 gp \n order 6: 12 gp']);
      end
      weight = parentArea*weight;        
    end

    function [coord, weight] = pointsTetra(nG)
      % Return Gauss integration points and weights for the reference tetrahedron
      % Reference tetrahedron: vertices at (0,0,0), (1,0,0), (0,1,0), (0,0,1)
      % Coordinates are in Cartesian coordinates
      % Weights are scaled by 1/6 (volume of reference tetrahedron)

      parentVolume = 1/6;
      switch nG
        case 1  % Degree 1 exact (1 point)
          coord = [1/4, 1/4, 1/4];   % Centroid
          weight = 1;

        case 4  % Degree 2 exact (4 points)
          a = 0.58541020;
          b = 0.13819660;
          coord = [
            b, b, b;
            a, b, b;
            b, a, b;
            b, b, a
            ];
          weight = repmat(1/4, 4, 1);

        case 5  % Degree 3 exact (5 points)
          coord = [
            1/4, 1/4, 1/4;
            1/2, 1/6, 1/6;
            1/6, 1/2, 1/6;
            1/6, 1/6, 1/2;
            1/6, 1/6, 1/6
            ];
          weight = [
            -4/5;
            9/20;
            9/20;
            9/20;
            9/20
            ];

        case 11  % Degree 4 exact (11 points) — Keast rule
          % 1 point at center
          coord1 = [1/4, 1/4, 1/4];
          w1 = (-74 / 5625);

          % 4 points at (a,a,a), (b,a,a), ..., permutations
          a = 11/14;
          b = 1/14;
          perm4 = [
            a, b, b;
            b, a, b;
            b, b, a;
            b, b, b
            ];
          w2 = (343 / 45000);

          % 6 points at (c,c,d), ..., permutations
          c = 0.399403576166799;
          d = 0.100596423833201;
          perm6 = [
            c, c, d;
            c, d, c;
            d, c, c;
            c, d, d;
            d, c, d;
            d, d, c
            ];
          w3 = (56 / 2250);

          coord = [
            coord1;
            perm4;
            perm6
            ];
          weight = [
            w1;
            repmat(w2, 4, 1);
            repmat(w3, 6, 1)
            ];

        otherwise
          error(['Unsupported number of gauss points for tetrahedron quadrature.' ...
            'Available integration order are: \n' ...
            ' order 1: 1 gp \n order 2: 4 gp \n order 3: 5 gp \n' ...
            ' order 4: 11 gp']);
      end
      weight = weight * parentVolume;        % refernce volume scaling
      assert(sum(weight)==parentVolume,'Wrong gauss point weight')
    end

  end
end

