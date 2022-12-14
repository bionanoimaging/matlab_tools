/* get the size for Khoros files 
   Bernd Rieger & Keith Lidke & Rainer Heintzmann
   June 2004-June 2006
*/

#include "mex.h"
#include <stdio.h>
#include <string.h>

void myexit(int arg) 
{ 
    mexErrMsgTxt("Fatal error in reading khoros file. Bailing out");
}

#include "khoros.h"


void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* format:
   writekhoros_info('blabla',dims,type)*/

/* nlhs # number left handed parameters
   plhs # left handed parameter array
   nrhs # number right handed parameters
   prhs # right handed parameter array*/
    ofstream from;
    char *TypeString; // will be filled in
    double *dp;
    char Var[2];
    char* matlabtype;
    char *input_buf;
    int   buflen,buflen2,status;

    //printf("int is %d\n",sizeof(int));
    //%printf("Long is %d\n",sizeof(long));
    //%printf("size_T is %d\n",sizeof(size_t));
    
    if (nrhs != 3) 
         mexErrMsgTxt("3 inputs required.");
    if (mxIsChar(prhs[0]) != 1)
         mexErrMsgTxt("Input must be a string.");
    if (mxGetM(prhs[0]) != 1)
      mexErrMsgTxt("Input must be a row vector.");
    if (mxIsChar(prhs[2]) != 1)
         mexErrMsgTxt("Input datatype must be a string.");
    if (mxGetM(prhs[2]) != 1)
      mexErrMsgTxt("Input datatype must be a row vector.");
    
    
    /* Get the length of the input string. */
    buflen = (mxGetM(prhs[0]) * mxGetN(prhs[0])) + 1;
    buflen2 = (mxGetM(prhs[2]) * mxGetN(prhs[2])) + 1;
    
    /* Allocate memory for input and output strings. */
    input_buf = (char*) mxCalloc(buflen, sizeof(char));
    matlabtype = (char*) mxCalloc(buflen2, sizeof(char));
    
    /* Copy the string data from prhs[0] into a C string */
    mxGetString(prhs[0], input_buf, buflen);
    mxGetString(prhs[2], matlabtype, buflen2);
    
    //printf("strlen %d,%d\n",strlen(input_buf),buflen);
    //printf("strlen %d, %d\n",strlen(matlabtype),buflen2);
    
    if (mxGetN(prhs[1])!=5)
      mexWarnMsgTxt("dims must be 5d vector");
    dp = mxGetPr(prhs[1]);
    
    if (strcmp(matlabtype,"uint8")==0) TypeString="Unsigned Byte";
    else if (strcmp(matlabtype,"sint32")==0) TypeString="Integer";
    else if (strcmp(matlabtype,"sfloat")==0) TypeString="Float";
    else if (strcmp(matlabtype,"dfloat")==0) TypeString="Double";
    else if (strcmp(matlabtype,"double")==0) TypeString="Double";
    else if (strcmp(matlabtype,"scomplex")==0) TypeString="Complex";
    else if (strcmp(matlabtype,"uint16")==0) TypeString="Unsigned Short";
    else if (strcmp(matlabtype,"dcomplex")==0) TypeString="Double Complex";   
    else if (strcmp(matlabtype,"bin")==0) TypeString="Unsigned Byte";   
    else {
      cerr << "Error: Datatype " << matlabtype << " unknown, cannot convert; bailing out!\n";
      return;
    }

    from.open(input_buf,ios::out | ios::binary);
    WriteKhorosHeader(&from,"Matlab2Khoros",TypeString, (int) dp[0], 
         (int) dp[1],(int)dp[2], (int) dp[3], (int) dp[4]);
    from.close();
    
    mxFree(input_buf);
    mxFree(matlabtype);
    return;
}
