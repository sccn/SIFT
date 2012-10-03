#include "mex.h"

/*
 * Arguments:
 * double *MatC       Pointer to matrix MatC, the
 * Kronecker product of the two matrix arguments in the
 * order they appear. The matrix C should be declared
 * double C[rowsa*rowsb][colsa*colsb] in the calling
 * routine.
 * double *A       Pointer to the first matrix in the Kronecker product.
 * int    rowsa    The number of rows of A.
 * int    colsa    The number of cols of A.
 * double *B       Pointer to the second matrix in the Kronecker product.
 * int    rowsb    The number of rows of B.
 * int    colsb    The number of cols of B.
 *
 * Return Values: Void                                                            //
 */


void kron(double *MatC, double *MatA, int rowsa, int colsa, double *MatB, int rowsb, int colsb) {
    int i, j, k, l, t, u, z=rowsa*rowsb, y=colsa*colsb;
    
    for (i = 0; i < z; i += rowsb){
        for (j = 0; j < y; j +=colsb){
            for (k = 0; k < rowsb; k++){
                for (l = 0; l < colsb; l++){
                    t = i/rowsb;
                    u = j/colsb;
                    *(MatC + (j +l)*z + i + k) =  *(MatB + l * rowsb + k) * *(MatA + u * rowsa + t); 
                }
            }
        }
    }
    
}


void mexFunction( int nlhs, mxArray *plhs[], int nrhs, const mxArray *prhs[]) {
    double *MatA, *MatB, *MatC;
    
    mwSize rowsa, colsa, rowsb, colsb;
    
    /* check for proper number of arguments */
    if(nrhs!=2)
        mexErrMsgTxt("Two inputs required.");
    else if(nlhs > 2)
        mexErrMsgTxt("Too many output arguments.");
    
    /*  create a pointer to the input vector r */
    MatA = mxGetPr(prhs[0]);
    MatB = mxGetPr(prhs[1]);
    
    /*  get the dimensions of the matrix input y */
    rowsa = mxGetM(prhs[0]);
    colsa = mxGetN(prhs[0]);
    rowsb = mxGetM(prhs[1]);
    colsb = mxGetN(prhs[1]);
    
    /*  set the output pointer to the output matrix */
    plhs[0] = mxCreateDoubleMatrix(rowsa * rowsb , colsa * colsb, mxREAL);
    
    /*  create a C pointer to a copy of the output matrix */
    MatC = mxGetPr(plhs[0]);
    
    /* call the C subroutine */
    kron(MatC, MatA, rowsa, colsa, MatB, rowsb, colsb);
    return;
}
