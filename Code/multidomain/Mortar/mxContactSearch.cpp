#include "mex.h"
#include <vector>

// Funzione ricorsiva per la traversal in tandem
void tandemTraversal(
    const double* treeNodes1, const double* treeNodes2,
    const mwSize* BBtree1, const mwSize* BBtree2,
    const mwSize* leaf2elem1, const mwSize* leaf2elem2,
    mwSize node1, mwSize node2,
    std::vector<mwSize>& rowIdx, std::vector<mwSize>& colIdx,
    mwSize nDirs)
{
    const double* node1box = &treeNodes1[node1 * 2 * nDirs];
    const double* node2box = &treeNodes2[node2 * 2 * nDirs];

    // Check intersection
    for (mwSize i = 0; i < nDirs; ++i) {
        double min1 = node1box[2*i];
        double max1 = node1box[2*i + 1];
        double min2 = node2box[2*i];
        double max2 = node2box[2*i + 1];
        if (min1 > max2 || min2 > max1) return;
    }

    bool isLeaf1 = (BBtree1[2*node1] == 0 && BBtree1[2*node1 + 1] == 0);
    bool isLeaf2 = (BBtree2[2*node2] == 0 && BBtree2[2*node2 + 1] == 0);

    if (isLeaf1 && isLeaf2) {
        rowIdx.push_back(leaf2elem1[node1] - 1); // C-style indexing
        colIdx.push_back(leaf2elem2[node2] - 1);
        return;
    }

    if (!isLeaf1) {
        tandemTraversal(treeNodes1, treeNodes2, BBtree1, BBtree2, leaf2elem1, leaf2elem2,
                        BBtree1[2*node1]-1, node2, rowIdx, colIdx, nDirs);
        tandemTraversal(treeNodes1, treeNodes2, BBtree1, BBtree2, leaf2elem1, leaf2elem2,
                        BBtree1[2*node1+1]-1, node2, rowIdx, colIdx, nDirs);
    } else {
        tandemTraversal(treeNodes1, treeNodes2, BBtree1, BBtree2, leaf2elem1, leaf2elem2,
                        node1, BBtree2[2*node2]-1, rowIdx, colIdx, nDirs);
        tandemTraversal(treeNodes1, treeNodes2, BBtree1, BBtree2, leaf2elem1, leaf2elem2,
                        node1, BBtree2[2*node2+1]-1, rowIdx, colIdx, nDirs);
    }
}

void mexFunction(int nlhs, mxArray* plhs[], int nrhs, const mxArray* prhs[])
{
    const mwSize* BBtree1 = mxGetUint64s(prhs[0]);
    const mwSize* BBtree2 = mxGetUint64s(prhs[1]);
    const double* treeNodes1 = mxGetDoubles(prhs[2]);
    const double* treeNodes2 = mxGetDoubles(prhs[3]);
    const mwSize* leaf2elem1 = mxGetUint64s(prhs[4]);
    const mwSize* leaf2elem2 = mxGetUint64s(prhs[5]);

    mwSize nDirs = mxGetN(prhs[2]) / 2;

    std::vector<mwSize> rowIdx;
    std::vector<mwSize> colIdx;

    tandemTraversal(treeNodes1, treeNodes2, BBtree1, BBtree2, leaf2elem1, leaf2elem2,
                    0, 0, rowIdx, colIdx, nDirs);

    // Output
    plhs[0] = mxCreateNumericMatrix(rowIdx.size(), 1, mxUINT64_CLASS, mxREAL);
    plhs[1] = mxCreateNumericMatrix(colIdx.size(), 1, mxUINT64_CLASS, mxREAL);
    mwSize* outRow = mxGetUint64s(plhs[0]);
    mwSize* outCol = mxGetUint64s(plhs[1]);
    for (mwSize i = 0; i < rowIdx.size(); ++i) {
        outRow[i] = rowIdx[i] + 1; // convert to MATLAB index
        outCol[i] = colIdx[i] + 1;
    }
}
