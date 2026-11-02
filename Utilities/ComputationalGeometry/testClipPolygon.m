function testClipPolygon()
% TESTCLIPPOLYGONMEX
% Simple graphical test for a polygon clipping MEX.
%
% Assumptions:
%   clip = clipPolygon(polySubject, polyClipper)
%
% Input polygons:
%   - Nx2 and Mx2 arrays
%   - vertices ordered consistently (CW or CCW)
%
% Output:
%   - Kx2 clipped polygon
%
% The script:
%   1) defines two test polygons
%   2) calls the clipping MEX
%   3) plots subject, clipper, and filled clipped region

clc;
close all;

% -------------------------------------------------------------------------
% Test polygons
% -------------------------------------------------------------------------
% Subject polygon
poly1 = [
    0.0 0.0
    3.0 0.5
    2.7 2.8
    1.2 3.4
   -0.2 2.0
];

% Clipper polygon
poly2 = [
    0.8 0.8
    3.2 1.0
    2.4 3.2
    0.6 2.6
];

% -------------------------------------------------------------------------
% Call clipping mex
% -------------------------------------------------------------------------
clip = clipPolygon(poly1, poly2);

% -------------------------------------------------------------------------
% Plot
% -------------------------------------------------------------------------
figure('Color','w');
hold on;
axis equal;
grid on;
box on;

title('Polygon clipping test');
xlabel('x');
ylabel('y');

% Plot subject polygon
plotClosedPolygon(poly1, 'k-', 1.5);
patch(poly1(:,1), poly1(:,2), [0.75 0.85 1.00], ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none');

% Plot clipper polygon
plotClosedPolygon(poly2, 'r-', 1.5);
patch(poly2(:,1), poly2(:,2), [1.00 0.80 0.80], ...
    'FaceAlpha', 0.20, 'EdgeColor', 'none');

% Plot clipped polygon
if ~isempty(clip) && size(clip,1) >= 3
    clip = orderPointsCCW(clip);
    patch(clip(:,1), clip(:,2), [0.20 0.70 0.30], ...
        'FaceAlpha', 0.65, 'EdgeColor', 'g', 'LineWidth', 2.0);
    plotClosedPolygon(clip, 'g-', 2.0);
else
    disp('No clipped polygon returned, or fewer than 3 vertices.');
end

legend({'Subject boundary', 'Subject fill', ...
        'Clipper boundary', 'Clipper fill', ...
        'Clipped fill', 'Clipped boundary'}, ...
        'Location','bestoutside');

% Vertex markers
plot(poly1(:,1), poly1(:,2), 'ko', 'MarkerFaceColor', 'k', 'MarkerSize', 5);
plot(poly2(:,1), poly2(:,2), 'ro', 'MarkerFaceColor', 'r', 'MarkerSize', 5);

if ~isempty(clip)
    plot(clip(:,1), clip(:,2), 'go', 'MarkerFaceColor', 'g', 'MarkerSize', 6);
end

hold off;

% -------------------------------------------------------------------------
% Print result
% -------------------------------------------------------------------------
disp('Subject polygon:');
disp(poly1);

disp('Clipper polygon:');
disp(poly2);

disp('Clipped polygon:');
disp(clip);

end


function plotClosedPolygon(P, lineSpec, lw)
% Plot polygon boundary closed
Q = [P; P(1,:)];
plot(Q(:,1), Q(:,2), lineSpec, 'LineWidth', lw);
end