function writeBCfiles_FixCenterPoints(cx, cy, cz, r0)
    % Fixing the particle center (ux=uy=uz=0)
    r = sqrt(cx.^2 + cy.^2 + cz.^2);
    constrainIndices = find(r < r0); % all nodes within that small sphere
    writeBCfiles('BCs/chemomech_u_0', 'NodeBC', 'Dir', {'Poromechanics', ...
        'x', 'y', 'z'}, 'Fixed center point', 0, 0, constrainIndices);
end