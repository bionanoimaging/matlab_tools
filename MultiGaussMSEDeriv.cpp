/* Compares a number of simulated Gaussian Spots with experimental data
 * as result, the error and the derivative with respect to all input vector components is returned
 * Examle : a=noise(51.2*exp(-((xx(20,20,10)-2).^2+(yy(20,20,10)-1.2).^2+(zz(20,20,10)+1.6).^2)/20)+33,'poisson')
           [fitted,params]=FitDataNDFast([0 1 0; 1 0 0],a)
 */

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// #define debug

// static int VaryOffset=1; 
// static int VarySlopeX=1;
// static int VarySlopeY=1;
// static int VaryInt=2; // 0: do not allow changes, 1: vary globally, 2: vary locally
// static int VaryPosX=2;
// static int VaryPosY=2;
// static int VaryPosZ=2;
// static int VaryWidthX=1;
// static int VaryWidthY=1;
// static int VaryWidthZ=1;

// bg sx sy px py wx wy int
#define numTotalParams(numDims,numSpots) (1+numDims+(numSpots)*(2*numDims+1))
#define numSpotsFromParams(numParams,numDims) (((double)(numParams-1-numDims)) / ((double)(2*numDims+1)))
#define PosIdx(aDim,gaussnum,numDims) (1+numDims+gaussnum*(2*numDims+1)+aDim)

// meaning of parameters below: 
// params: first global parameters, then blocks of size paramSize. One block for each particle
// Function is:   offset + sx*x + sy*y + brightness * exp(-( (x-x0)/wx)^2 + ...)
// [0]: global offset
// [1 .. numDims]: global slopes
// [numDims+1 .. 2*numDims+1]: position along the dimension
// [2*numDims+1 .. 3*numDims+1]: width along the dimension
// [3*numDims+2]: brightness (pre-exponential factor)
double dosim(double * params, int numDims, int numSpots, int * pos, double * resDeriv)   // Calculates a Gaussian and its derivatives
{
    int paramIdx=0,d,n;
    double result=params[paramIdx++],sq,ssq,myExp,tmp;
    resDeriv[0]=1;  // derivative of global offset
    for (d=0;d<numDims;d++)  // global slope
    {
        resDeriv[paramIdx]=(double) pos[d];  // derivative of global slopes
        result += pos[d]*params[paramIdx];
        paramIdx++;
    }
    //printf("PosX : %d, y %d\n",pos[0],pos[1]);
    for (n=0;n<numSpots;n++)
    {
        ssq=0;
        for (d=0;d<numDims;d++)  // position parameters
        {sq=(pos[d]-params[paramIdx+d]);  // accesses the position parameters
         sq /= params[paramIdx+numDims+d];  // accesses the width parameters
        //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
         ssq +=sq*sq;}
        myExp=exp(-ssq);  // derivative of global slopes
        for (d=0;d<numDims;d++)  // To fill in the derivatives
        {   tmp = 2*params[paramIdx+2*numDims]*myExp*(pos[d]-params[paramIdx+d])/params[paramIdx+numDims+d]/params[paramIdx+numDims+d];
            resDeriv[paramIdx+d] = tmp; // position derivatives
            resDeriv[paramIdx+numDims+d] = tmp*(pos[d]-params[paramIdx+d])/params[paramIdx+numDims+d];  // width derivatives
        }
        paramIdx+=2*numDims;
        resDeriv[paramIdx]=myExp;  // intensity derivative
        result += params[paramIdx]*myExp;  // forward simulation
        paramIdx++;
        //printf("Intensity : %g\n",params[paramSize*n]);
    }
    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    return result;
}

// this function computes the meas square error comparing data with simulation
// params is the array of spot positions
// paramsg are the global parameters
double do_mse(double * mydata, int * sizes, double * params, int numDims, int numSpots, 
        double * res, // result of simulation
        double * resdiff,  // difference to data
        double * resDeriv,  // tmp vector to compute derivative
        double * resgrad)  // gradient vector (all coordinates)
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    //for (int y=0;y<numSpots;y++) 
    //   for (int x=0;x<paramSpots;x++) 
    //      printf("MSE: Param %d, %d : %g\n",x,y,params[x+y*paramsize]);
    double result=0,tmp;
    int i=0,pos[3],d;
//    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
//        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
//            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
    for (pos[2]=-floor((double) sizes[2]/2);pos[2]< -floor((double) sizes[2]/2)+sizes[2];pos[2]++)
    for (pos[1]=-floor((double) sizes[1]/2);pos[1]< -floor((double) sizes[1]/2)+sizes[1];pos[1]++)
    for (pos[0]=-floor((double) sizes[0]/2);pos[0]< -floor((double) sizes[0]/2)+sizes[0];pos[0]++)
                {
                tmp = dosim(params, numDims, numSpots, pos, resDeriv);                 
                if (res != 0)
                    res[i] = tmp; // save simulation
                tmp=mydata[i]-tmp; // residuum
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference
                for (d=0;d<numTotalParams(numDims,numSpots);d++)
                    resgrad[d] += - 2 * tmp * resDeriv[d];  // gradient of mse
                result += tmp*tmp;  // summing error terms
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

double do_idiv(double * mydata, int * sizes, double * params, int numDims, int numSpots, 
        double * res, // result of simulation
        double * resdiff,  // residual (here ratio)
        double * resDeriv,  // tmp vector to compute derivative
        double * resgrad)  // gradient vector
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp;
    int i=0,pos[3],d;
//     for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
//         for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
//             for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
    for (pos[2]=-floor((double) sizes[2]/2);pos[2]< -floor((double) sizes[2]/2)+sizes[2];pos[2]++)
    for (pos[1]=-floor((double) sizes[1]/2);pos[1]< -floor((double) sizes[1]/2)+sizes[1];pos[1]++)
    for (pos[0]=-floor((double) sizes[0]/2);pos[0]< -floor((double) sizes[0]/2)+sizes[0];pos[0]++)
                {
                tmp = dosim(params, numDims, numSpots, pos, resDeriv);                 
                if (res != 0)
                    res[i] = tmp; // save idiv image
                for (d=0;d<numTotalParams(numDims,numSpots);d++)
                    resgrad[d] += resDeriv[d]*(1.0 - mydata[i] / tmp);
                if (mydata[i] !=0)
                    tmp=mydata[i]*log(mydata[i]/tmp)-(mydata[i]-tmp);  // i-divergence with sterling's approximation
                else
                    tmp=-(mydata[i]-tmp);
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

double do_fidiv(double * mydata, int * sizes, double * params, int numDims, int numSpots, 
        double * res, 
        double * resdiff,
        double * resDeriv,  // tmp vector to compute derivative
        double * resgrad)  
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp;
    int i=0,pos[3],d;
//     for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
//         for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
//             for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
    for (pos[2]=-floor((double) sizes[2]/2);pos[2]< -floor((double) sizes[2]/2)+sizes[2];pos[2]++)
    for (pos[1]=-floor((double) sizes[1]/2);pos[1]< -floor((double) sizes[1]/2)+sizes[1];pos[1]++)
    for (pos[0]=-floor((double) sizes[0]/2);pos[0]< -floor((double) sizes[0]/2)+sizes[0];pos[0]++)
                {
                tmp = dosim(params, numDims, numSpots, pos, resDeriv);                 
                if (res != 0)
                    res[i] = tmp; // save idiv image
                for (d=0;d<numTotalParams(numDims,numSpots);d++)
                    resgrad[d] += resDeriv[d]*(1.0 - mydata[i] / tmp);
                tmp=tmp-mydata[i]*log(tmp);  // fast version omitting constants
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    return result;
}

void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
/* format:
   writekhoros_info('blabla',dims,type)
   MultiGaussMSE(Data)  to prepare with experimental data
   MultiGaussMSE(Matrix with parameters, Global parameters) each row containing intensity and positions */

/* nlhs # number left handed parameters
   plhs # left handed parameter array
   nrhs # number right handed parameters
   prhs # right handed parameter array*/
    double *data=0;  // Let's hope this one survives the individual calls
    double * params, * fixedparams, result, * res=0, * resdiff=0;
    char * input_buf;
    int PsizeX,numSpots, buflen; // GPsize, sizes[5],
    static double * mydata=0;   // stores the data to fit
//    static double * cparams=0;  // 
    static double * resDeriv=0;  // stores temporarily the derivatives of the Gaussian with respect to all the parameters 
    static double * resgrad=0;  // stores the gradient vector (which is then transferred to Matlab)
    static int allocated=0;
    static int sizes[100],nd=1,totalsize=1,numdims;   // To hell with ppl who use more than 100 dimensions!
    enum method {mse, idiv,fidiv};
    static enum method mymethod=mse;
    
    if (nrhs != 1 && nrhs != 3) 
    {
        printf ("Nr of parameters: %d\n",nrhs);
         mexErrMsgTxt("1 or 3 inputs required");
    }
    //if (mxIsChar(prhs[0]) != 1)
    //     mexErrMsgTxt("Input must be a string.");


    if (nrhs > 1)  // archieve the data by saving the pointer and remember which method to use
    {
        int ii,i;
        const int *sz;
    if (nrhs < 3)  
            mexErrMsgTxt("When submitting data, three arguments are required: data, method-string, dimensions!");
        
    /* Get the length of the input string. */
        sz = mxGetDimensions(prhs[0]);
        nd = mxGetNumberOfDimensions(prhs[0]);
        totalsize=1;
#ifdef debug
        printf("Data dims %d\n",nd);
#endif
        for (ii=0;ii<100;ii++)
        {
            sizes[ii]=1;
        }

        for (ii=0;ii<nd;ii++)
        {
           sizes[ii]=sz[ii];  // save these values
#ifdef debug
            printf("dim %d, size %d\n",ii,sizes[ii]);
#endif
           totalsize *= sz[ii];
        }
        if (totalsize == 1)
            mexErrMsgTxt("Data array contains only 1 number. Probably the data was submitted as DipImage but needs to be converted to double!");
            
        data = mxGetPr(prhs[0]);
        if (mydata != 0) free(mydata);
        mydata=(double *) calloc(totalsize,sizeof(double));
        for (i=0;i<totalsize;i++) 
            mydata[i]=data[i];
        if (nlhs != 0)
            mexErrMsgTxt("When submitting data, no output is returned!");
        // printf("Bufflen is %d x %d, pointer is %x, copied to %x\n",dataSizeX,dataSizeY,data,mydata);
    /* Get the length of the input string. */
    buflen = (int) (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 4;
    
    /* Allocate memory for input and output strings. */
    input_buf = (char*) mxCalloc(buflen, sizeof(char));
    
    /* Copy the string data from prhs[0] into a C string */
    mxGetString(prhs[1], input_buf, buflen);

    if (strcmp(input_buf,"mse") == 0)
        mymethod = mse;
    else if (strcmp(input_buf,"idiv") == 0)
        mymethod = idiv;
    else if (strcmp(input_buf,"fidiv") == 0)
        mymethod = fidiv;
    else
    {
        printf("Requested method was %s\n",input_buf);
        mexErrMsgTxt("Invalid method. Valid methods are : 'mse', 'idiv' and 'fidiv'");
    }
    /* Copy the string data from prhs[0] into a C string */
    numdims=(int) (* mxGetPr(prhs[2]));

    mxFree(input_buf);
    }
    else   // Only one right side argument was given -> argument vector
    {
        PsizeX = (int) mxGetM(prhs[0]);
        if (mxGetN(prhs[0]) > 1)
            mexErrMsgTxt("All parameters should be in a single vector.(try using the transpose)");
        if (numSpotsFromParams(PsizeX,numdims) != (double) ((int) numSpotsFromParams(PsizeX,numdims)))
            mexErrMsgTxt("Number of parameters (globals and rest) does not match with number of dimensions to fit. Parameters are Offset, SlopeX, SlopeY,.., {PosX,PosY,..,WidthX,WidhtY, ..., Int}");
        else
        {
            numSpots = (int) numSpotsFromParams(PsizeX,numdims);
        }
        // printf("PsizeX %d, numSpots %g, %g\n",PsizeX,numSpotsFromParams(PsizeX,numdims),(double) ((int) numSpotsFromParams(PsizeX,numdims)) );
        params = mxGetPr(prhs[0]);

        if (allocated < numSpots)
        {
            //if (cparams != 0)
            //    free(cparams);
            //cparams=(double *) calloc(numTotalParams(numdims,numSpots),sizeof(double));
            if (resDeriv != 0)
                free(resDeriv);
            resDeriv=(double *)calloc(numTotalParams(numdims,numSpots),sizeof(double));
            if (resgrad != 0)
                free(resgrad);
            resgrad=(double *) calloc(numTotalParams(numdims,numSpots),sizeof(double));
            allocated=numTotalParams(numdims,numSpots);
#ifdef debug
           printf("Allocated arrays, NumTotalParams: %d\n",numTotalParams(numdims,numSpots));
#endif
        }

        //if (mxGetM(prhs[1]) != 1)
        //    mexErrMsgTxt("Global parameters must be a row vector.");
    

        //GPsize = mxGetN(prhs[1]);
        //fixedparams = mxGetPr(prhs[1]);
        //printf("ParamsizeX: %d, Y %d\n",PsizeX,numSpots);
        //printf("dataSizeX: %d, Y %d\n",dataSizeX,dataSizeY);
        if (mydata != 0) {
            //for (int i=0;i<dataSizeX*dataSizeY;i++) 
            //    printf("%d: %g\n",i,mydata[i]);
            //for (int y=0;y<numSpots;y++) 
            //  for (int x=0;x<PsizeX;x++) 
            //    printf("Param %d, %d : %g\n",x,y,params[x+y*PsizeX]);

            
//            for (int d=0;d<numTotalParams(numdims,numSpots);d++)   // copy params to cparams
//                    cparams[d] = params[d];  

//             for (int d=0;d<numSpots;d++)   // to account for center corresponding to zero
//                 for (int n=0;n<numdims;n++)   // to account for center corresponding to zero
//                 {
//                     cparams[PosIdx(n,d,numdims)] = params[PosIdx(n,d,numdims)] + floor((double) sizes[n]/2);  
//                     printf("Changing Spot %d,dim %d, from %g to %g\n",d,n,params[PosIdx(n,d,numdims)],cparams[PosIdx(n,d,numdims)] );
//                 }
            int d;
            if (nlhs >= 3)
            {
                plhs[2] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                res = mxGetPr(plhs[2]);  // result array                
            }
            else res=0;
            
            if (nlhs >= 4)
            {
                plhs[3] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                resdiff = mxGetPr(plhs[3]);  // residual aray
            }
            else resdiff=0;
            //for (int y=0;y<numSpots;y++) 
            //  for (int x=0;x<PsizeX;x++) 
            //    printf("Param %d, %d : %g\n",x,y,params[x+y*PsizeX]);
#ifdef debug
           printf("NumSpots %d, numDims %d, nd %d, NumParams %d\n",numSpots, numdims, nd,numTotalParams(numdims,numSpots));
#endif
            for (d=0;d<numTotalParams(numdims,numSpots);d++)   // clear the gradient as this is computed as a sum
                resgrad[d]=0;
            
            if (PsizeX != numTotalParams(numdims,numSpots))
                    mexErrMsgTxt("number of dimension does not match vector length");
            
            switch (mymethod)
            {
                case mse:
                    result=do_mse(mydata, sizes, params, numdims, numSpots, res, resdiff, resDeriv, resgrad); 
                break;
                case idiv:
                    result=do_idiv(mydata, sizes, params, numdims, numSpots, res, resdiff, resDeriv, resgrad); 
                break;
                case fidiv:
                    result=do_fidiv(mydata, sizes, params, numdims, numSpots, res, resdiff, resDeriv, resgrad); 
                break;
                default:
                    mexErrMsgTxt("Undefined method. Valid methods are : 'mse', 'idiv' and 'fidiv'");
            }
            if (nlhs >= 1)
            {
                double * dp;
                plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
                dp = mxGetPr(plhs[0]);
                (* dp) = result;
            }
            if (nlhs >= 2)
            {
                const int *sz;
                double * dp;
                int nd = mxGetNumberOfDimensions(prhs[0]),d;
                sz = mxGetDimensions(prhs[0]);
                plhs[1] = mxCreateNumericArray(nd, sz, mxDOUBLE_CLASS, mxREAL);
                // plhs[1] = mxCreateDoubleMatrix(numTotalParams(numdims,numSpots),1, mxREAL);
                dp = mxGetPr(plhs[1]);
#ifdef debug                
                printf("NumTotalParams: %d\n",numTotalParams(numdims,numSpots));
#endif
                for (d=0;d<numTotalParams(numdims,numSpots);d++)   // transfer the gradient to Matlab
                    dp[d] = resgrad[d];
            }
         // mexErrMsgTxt("Point 3");return;
#ifdef debug
         printf("Result %g\n",result);
#endif
        }
        else
            mexErrMsgTxt("Please provide just the data matrix first");
    }
    
    return;
}
