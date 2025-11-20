function writeBCfiles_PointwiseConstraint(cx, cy, cz, tol)
    % Fix the displacement at the center (u = [0 0 0])
    Point_000 = find(all(abs([cx-0, cy-0 cz-0]) < tol,2));
    writeBCfiles('BCs/chemomech_u_0', 'NodeBC', 'Dir', {'Poromechanics', ...
        'x', 'y', 'z'}, 'Fix_center', 0, 0, Point_000);
    
    % Fix ux = 0 at [0 0 1]
    Point_001 = find(all(abs([cx-0, cy-0 cz-1]) < tol,2));
    writeBCfiles('BCs/chemomech_u_x', 'NodeBC', 'Dir', {'Poromechanics', ...
        'x'}, 'Fix_ux', 0, 0, Point_001);

    % % Fix ux = 0 at [0 1 0]
    % Point_010 = find(all(abs([cx-0, cy-1 cz-0]) < tol,2));
    % writeBCfiles('BCs/chemomech_u_x', 'NodeBC', 'Dir', {'Poromechanics', ...
    %     'x'}, 'Fix_ux', 0, 0, Point_010);
    
    % Fix uy=uz=0 at [1 0 0]
    Point_100 = find(all(abs([cx-1, cy-0 cz-0]) < tol,2));
    writeBCfiles('BCs/chemomech_u_yz', 'NodeBC', 'Dir', {'Poromechanics', ...
        'y', 'z'}, 'Fix_uy_uz', 0, 0, Point_100);
end