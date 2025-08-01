%% Pore Network Generation and Analysis
% This script simulates a 3D porous medium and performs analysis on the generated pore structure.
% License: MIT License
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, including without limitation the rights
% to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
% copies of the Software, and to permit persons to whom the Software is
% furnished to do so, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
% 
% THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
% IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
% FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
% AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
% LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
% OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
% SOFTWARE.
%
% Contact Information:
% Author: Jimmy X. Li
% Email: Jimmy.Li@uq.edu.au
% Affiliation: The University of Queensland
%
% Citation:
% If you use this code in your research, please cite the following :
% Jimmy Li (2024). Pore Network Generation and Analysis 
% (https://www.mathworks.com/matlabcentral/fileexchange/170286-pore-network-generation-and-analysis), 
% MATLAB Central File Exchange. Retrieved July 26, 2024.
close all
clc
clear all
tic
% Define the 3D grid size
nx = 100;
ny = 100;
nz = 100;
% Initialize the 3D matrix to represent the domain (1-pore, 0-solid)
domain = zeros(nx, ny, nz);
% Target porosity and average pore size
targetPorosity = 0.75; % Adjust as needed
targetPoreSize = 20; % Target average pore size (e.g., radius)
% Variables to track total pore volume and number of pores
totalPoreVolume = 0;
numberOfPores = 0;
% Seed for random number generator for reproducibility
rng(0);
while sum(domain(:)) / numel(domain) < targetPorosity
    % Calculate current average pore size
    currentAveragePoreSize = 0;
    if numberOfPores > 0
        currentAveragePoreSize = (3/4/pi) * (totalPoreVolume / numberOfPores)^(1/3);
    end
    
    % Adjust growth size based on the difference from the target pore size
    growthSize = max(1, round(targetPoreSize - currentAveragePoreSize + 1));
    
    % Select a random starting point
    seedX = randi([1 nx]);
    seedY = randi([1 ny]);
    seedZ = randi([1 nz]);
    
    if domain(seedX, seedY, seedZ) == 0
        % Initialize temporary storage for the new pore
        newPoreVolume = 0;
        
        for x = max(1, seedX-growthSize):min(nx, seedX+growthSize)
            for y = max(1, seedY-growthSize):min(ny, seedY+growthSize)
                for z = max(1, seedZ-growthSize):min(nz, seedZ+growthSize)
                    if sqrt((x-seedX)^2 + (y-seedY)^2 + (z-seedZ)^2) <= growthSize
                        if domain(x, y, z) == 0 % Only if it's solid
                            domain(x, y, z) = 1; % Assign as pore space
                            newPoreVolume = newPoreVolume + 1;
                        end
                    end
                end
            end
        end
        
        % Update the total pore volume and number of pores
        if newPoreVolume > 0
            totalPoreVolume = totalPoreVolume + newPoreVolume;
            numberOfPores = numberOfPores + 1;
        end
    end
end
% Example usage
% Assume 'domain' is your 3D binary matrix
% domain = rand(100, 100, 100) > 0.75; % Example: generate a random domain
[numCells, numFaces, numEdges, numVertices, vertices] = analyzeDomain(domain);
fprintf('Number of Cells (Voxels): %d\n', numCells);
fprintf('Number of Faces: %d\n', numFaces);
fprintf('Number of Edges: %d\n', numEdges);
fprintf('Number of Vertices: %d\n', numVertices);
% Vertices are stored in the 'vertices' variable
% Visualization (not repeated for brevity; see previous example)
% Assuming 'domain' contains your 3D porous medium data where 1 represents pores
% Find the isosurface for the boundary between solid and pore spaces
[f,v] = isosurface(domain, 0.5); % Adjust the isovalue if needed
[f_,v_] = isosurface(1-domain, 0.5); % Adjust the isovalue if needed
% Create a figure to visualize the mesh
% Visualization of the skeleton can be complex due to the 3D nature.
% Here's a basic approach using volumeViewer or sliceViewer for an interactive exploration
volumeViewer(domain);
volumeViewer(1-domain);
% or use sliceViewer for a slice-by-slice exploration
% sliceViewer(poreSkeleton);
% Plot the skeleton points
figure (1)
%pores
patch('Faces',f,'Vertices',v,'FaceColor','k','EdgeColor','none');
% Improve the lighting and appearance
camlight right; lighting phong; axis equal; 
hold on 
%solid
patch('Faces',f_,'Vertices',v_,'FaceColor','r','EdgeColor','none');
axis equal;
view(3)
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Mesh of Solid Part in Porous Medium');
figure(2)
% Assuming 'domain' is your 3D binary image with 1s for pores and 0s for solid
% Perform skeletonization
poreSkeleton = bwskel(logical(domain));
% Find the coordinates of the skeleton points
[skeletonX, skeletonY, skeletonZ] = ind2sub(size(poreSkeleton), find(poreSkeleton));
patch('Faces',f,'Vertices',v,'FaceColor','blue','EdgeColor','none');
% Improve the lighting and appearance
camlight right; lighting phong; axis equal; 
hold on 
plot3(skeletonX, skeletonY, skeletonZ, 'r.');
axis equal;
view(3)
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Mesh of Solid Part in Porous Medium');
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Assuming 'domain' is your 3D binary image with 1s for pores
EDT = bwdist(~domain);
% Flatten the EDT array to a vector for histogram analysis
EDT_vector = EDT(:);
% Filter out zero distances (i.e., solid matrix points)
EDT_vector = EDT_vector(EDT_vector > 0);
% Convert distances to diameters by multiplying by 2
poreDiameters = EDT_vector * 2;
% Generate histogram data
% Adjust the number of bins or specify bin edges for finer control over the plot
[n, edges] = histcounts(log10(poreDiameters), 50);
% Plot pore size histogram on a log scale
figure (3)
semilogx(edges(1:end-1), n, 'b*-');
xlabel('Log10(Pore Diameter)');
ylabel('Frequency');
title('Pore Size Distribution on Log Scale');
grid on;
%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% Exportation
% Assuming 'domain' is your 3D binary matrix
% Generate the isosurface; adjust the isovalue as needed
[F,V] = isosurface(domain, 0.5); % F: faces, V: vertices
% Call stlwrite with the filename and isosurface data
stlwrite('porous_media.stl', F, V);
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
cylindricalVolume = cylinderify(domain);
rcylindricalVolume = cylinderify(1-domain);
% Find the isosurface for the boundary between solid and pore spaces
[f_cylinder,v_cylinder] = isosurface(cylindricalVolume, 0.5); % Adjust the isovalue if needed
[f_cylinder_r,v_cylinder_r] = isosurface(rcylindricalVolume, 0.5); % Adjust the isovalue if needed
% Create a figure to visualize the mesh
figure (4)
patch('Faces',f_cylinder,'Vertices',v_cylinder,'FaceColor','k','EdgeColor','none');
hold on
patch('Faces',f_cylinder_r,'Vertices',v_cylinder_r,'FaceColor','y','EdgeColor','r');
% Improve the lighting and appearance
camlight right; lighting phong; axis equal; 
xlabel('X'); ylabel('Y'); zlabel('Z');
title('3D Mesh of Solid Part in Porous Medium');
view(3)
volumeViewer(cylindricalVolume);
volumeViewer(rcylindricalVolume);
% Assuming 'domain' is your binary 3D matrix where fractures are indicated by 1s
totalSurfaceArea = calculateFractureSurfaceArea(domain);
disp(['Total surface area of fractures: ', num2str(totalSurfaceArea)]);
toc
function cylindricalVolume = cylinderify(domain)
% export cylinderical shape
% Domain dimensions
[nx, ny, nz] = size(domain);
% Cylinder specifications
cylinderRadius = min(nx, ny) / 4; % Example radius
cylinderHeight = nz; % Example height, spanning the entire Z dimension
% Cylinder center (assuming it's centered in the XY plane)
centerX = nx / 2;
centerY = ny / 2;
% Initialize a mask for the cylindrical volume (0s and 1s, 1s inside the cylinder)
cylinderMask = zeros(size(domain));
% Loop through each point in the domain to determine if it's inside the cylinder
for x = 1:nx
    for y = 1:ny
        for z = 1:nz
            % Calculate the distance from the point to the cylinder's axis
            distanceFromCenter = sqrt((x - centerX)^2 + (y - centerY)^2);
            % Check if the point is within the cylinder's radius and height
            if distanceFromCenter <= cylinderRadius && z <= cylinderHeight
                cylinderMask(x, y, z) = 1;
            end
        end
    end
end
% Apply the mask to select the cylindrical volume from the domain
cylindricalVolume = domain .* cylinderMask;
end
function totalSurfaceArea = calculateFractureSurfaceArea(domain)
    % Calculate the surface area of the fractures within the 3D domain.
    % The domain is a binary 3D matrix where 1 represents fracture space.
    % Pad the domain with zeros on all sides to simplify boundary checks
    paddedDomain = padarray(domain, [1 1 1], 'both');
    % Initialize surface area count
    totalSurfaceArea = 0;
    % Voxel dimensions - assuming cubic voxels for simplicity, with a unit length of 1
    voxelFaceArea = 1; % Adjust this if your voxels are not unit cubes
    % Iterate through each voxel in the padded domain
    [nx, ny, nz] = size(paddedDomain);
    for x = 2:nx-1
        for y = 2:ny-1
            for z = 2:nz-1
                if paddedDomain(x, y, z) == 1  % Check only fracture voxels
                    % Check all six neighbors (face sharing)
                    if paddedDomain(x+1, y, z) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                    if paddedDomain(x-1, y, z) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                    if paddedDomain(x, y+1, z) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                    if paddedDomain(x, y-1, z) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                    if paddedDomain(x, y, z+1) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                    if paddedDomain(x, y, z-1) == 0
                        totalSurfaceArea = totalSurfaceArea + voxelFaceArea;
                    end
                end
            end
        end
    end
end
function stlwrite(filename, Faces, Vertices)
    % A simple function to write STL files
    % Inputs:
    % - filename: String, name of the file to create
    % - Faces: m x 3 matrix of indices that form triangular faces
    % - Vertices: n x 3 matrix of vertex coordinates
    % Open the file for writing
    fid = fopen(filename, 'w');
    
    % Write the STL header
    fwrite(fid, zeros(1, 80), 'uint8'); % 80-byte header
    fwrite(fid, size(Faces, 1), 'uint32'); % Number of triangles
    
    % Write the face data
    for i = 1:size(Faces, 1)
        % Normal vector (not used but required by the format)
        fwrite(fid, [0 0 0], 'float32');
        
        % Vertex coordinates of the triangle
        fwrite(fid, Vertices(Faces(i,:), :)', 'float32');
        
        % Attribute byte count
        fwrite(fid, 0, 'uint16');
    end
    
    % Close the file
    fclose(fid);
end
function [numCells, numFaces, numEdges, numVertices, vertices] = analyzeDomain(domain)
    % Extract the isosurface from the binary domain
    [faces, vertices] = isosurface(domain, 0.5);
    
    % Number of vertices is directly obtained
    numVertices = size(vertices, 1);
    
    % Number of faces is directly obtained from the isosurface
    numFaces = size(faces, 1);
    
    % Edges need to be calculated from the face-vertex connectivity
    edgeMap = containers.Map('KeyType', 'char', 'ValueType', 'int32');
    numEdges = 0;
    for i = 1:numFaces
        face = faces(i, :);
        % Create a list of edges for this face, ensuring uniqueness
        edges = nchoosek(face, 2);
        for j = 1:size(edges, 1)
            edgeStr = sprintf('%d_%d', min(edges(j,:)), max(edges(j,:)));
            if ~isKey(edgeMap, edgeStr)
                edgeMap(edgeStr) = 1;
                numEdges = numEdges + 1;
            end
        end
    end
    
    % For a binary domain represented as a 3D grid, cells are the grid elements (voxels)
    % Here, cells are considered as the non-zero elements in the domain
    numCells = nnz(domain);  % nnz counts the number of non-zero elements
    
    % Output the vertices as well
end
