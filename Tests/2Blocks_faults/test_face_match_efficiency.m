function test_face_match_efficiency()
    % Define face-node map for hexahedra
    faceNodes = [
        1 2 3 4;
        1 4 5 8;
        1 2 5 6;
        2 3 6 7;
        3 4 7 8;
        5 6 7 8
    ];

    sortedFaceNodes = sort(faceNodes, 2);
    nTests = 5e5;
    rng(0);  % For reproducibility

    % Generate random permutations of faceNodes for testing
    testFaces = zeros(nTests, 4);
    faceIndices = randi(6, nTests, 1);  % Random valid face indices

    for i = 1:nTests
        testFaces(i, :) = faceNodes(faceIndices(i), randperm(4));  % Random permutation
    end

    % ----------- Method 1: ismember with 'rows' -----------
    tic
    for i = 1:nTests
        sortedTest = sort(testFaces(i, :));
        id1 = find(ismember(sortedFaceNodes, sortedTest, 'rows'), 1);
    end
    t1 = toc;

    % ----------- Method 2: bsxfun/all -----------
    tic
    for i = 1:nTests
        sortedTest = sort(testFaces(i, :));
        id2 = find(all(bsxfun(@eq, sortedFaceNodes, sortedTest), 2), 1);
    end
    t2 = toc;

    fprintf("ismember time: %.4f s\n", t1);
    fprintf("bsxfun time:   %.4f s\n", t2);
    fprintf("Speed-up:      %.2fx\n", t1 / t2);
end