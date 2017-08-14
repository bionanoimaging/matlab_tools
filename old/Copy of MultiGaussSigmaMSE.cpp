/* Compares a number of simulated Gaussian Spots with experimental data

 * Examle : a=exp(-((xx(20,20)-2).^2+(yy(20,20)-1.2).^2)/20)
           [fitted,params]=FitDataNDFast([0 1 0; 1 0 0],a)
 */

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

// The function below calculates the N-dimensional Gaussian value at a given position pos
double dosim(double * params, int paramSize, int numparams, double * gparams, int * pos)
{
    double result=gparams[0],sq,ssq;  // offset is accounted for first gparams[0]
    //printf("PosX : %d, y %d\n",pos[0],pos[1]);
    for (int n=0;n<numparams;n++)
    {
        ssq=0;
        int numdim=(paramSize-1)/2;
        for (int d=0;d<numdim;d++)
        {sq=(pos[d]-params[paramSize*n+d+1])/params[paramSize*n+numdim+d+1];  // positions are in [paramSize*n+d+1], sigmas in [paramSize*n+numdim+d+1]
        //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
         ssq +=sq*sq;}
        result += params[paramSize*n]*exp(-ssq);                              // intensities are in params[paramSize*n]
        //printf("Intensity : %g\n",params[paramSize*n]);
    }
    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    return result;
}

// The function below calculates the Value (return value) and the Gradient of a pixel at position pos
double doGradSim(double * mygradient, double * params, int paramSize, int numparams, double * gparams, int * pos)
{
    double result=0.0,spotresult,sq,ssq;  
    //printf("PosX : %d, y %d\n",pos[0],pos[1]);
    for (int n=0;n<numparams;n++)  // goes through all the spots
    {
        ssq=0;
        int numdim=(paramSize-1)/2;   // dimension dependent are intensities and sigmas
        for (int d=0;d<numdim;d++)
        {sq=(pos[d]-params[paramSize*n+d+1])/params[paramSize*n+numdim+d+1];  // positions are in [paramSize*n+d+1], sigmas in [paramSize*n+numdim+d+1]
         //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
         ssq +=sq*sq;}
        spotresult = exp(-ssq);                              // intensities are in params[paramSize*n]
        mygradient[1+paramSize*n] = spotresult;                // gradient with respect to I
        spotresult *= params[paramSize*n];                   // mulitply by I
        result += spotresult;
        for (int d=0;d<numdim;d++)
        { sq= 2.0*(pos[d]-params[paramSize*n+d+1])/(params[paramSize*n+numdim+d+1] * params[paramSize*n+numdim+d+1]);
          mygradient[1+paramSize*n+d+1] = spotresult * sq ;   // positional gradient
          mygradient[1+paramSize*n+numdim+d+1] = spotresult * sq * (pos[d]-params[paramSize*n+d+1])/params[paramSize*n+numdim+d+1];   // sigma gradient
        }
        //printf("Intensity : %g\n",params[paramSize*n]);
    }
         result += gparams[0];               // offset is accounted for first gparams[0]
         // gradient with respect to global background is always 1.0
         mygradient[0] = 1.0;
    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    return result;
}

inline bool nextPos(dimensions,pos)
{
    if (mask)
    {
        pixel++;  //global pixel index
        pos=&mask[pixel*dimensions]; // mask stores all the X, Y Z values that are in the mask
    }
    else  // use ordinary indexing
    {
        d=0;
        pixel++;
        while ((++pos[d])>=sizes[d])
        {
            pos[d++]=0;
            if (d>= dimensions)
                return false;
        }
    }
    return true;
}

// this function computes the meas square error comparing data with simulation
// params is the array of spot positions
// paramsg are the global parameters
double do_mse(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double *resdiff, double *grad)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    //for (int y=0;y<numparams;y++) 
    //   for (int x=0;x<paramsize;x++) 
    //      printf("MSE: Param %d, %d : %g\n",x,y,params[x+y*paramsize]);
    double result=0,tmp,tmp2,* tmpgrad=0;
    int i=0,pos[3];
    if (grad != 0)
    {
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (int s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                if (grad != 0) {
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    tmp2=mydata[i]-tmp;
                    for (int s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (-2.0)*tmpgrad[s]*tmp2;  // sum over all the pixels
                }
                else
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 
                if (res != 0)
                    res[i] = tmp; // save simulation
                tmp=mydata[i]-tmp;
                if (resdiff != 0)
                    resdiff[i] = tmp; // save difference                
                result += tmp*tmp;     // sum of squares over all pixels
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    if (tmpgrad != 0)
        free(tmpgrad);
    return result;
}

double do_idiv(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff, double *grad)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp,*tmpgrad=0;
    int i=0,pos[3];
    if (grad != 0)
    {
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (int s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                if (grad != 0) {
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    for (int s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (1.0-mydata[i]/tmp)*tmpgrad[s];  // sum over all the pixels
                }
                else
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 

                if (res != 0)
                    res[i] = tmp; // save idiv image
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
    if (tmpgrad != 0)
        free(tmpgrad);
    return result;
}

double do_fidiv(double * mydata, int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff, double *grad)
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp,*tmpgrad=0;
    int i=0,pos[3];
    if (grad != 0)
    {
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (int s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    for (pos[2]=0;pos[2]< sizes[2];pos[2]++)
        for (pos[1]=0;pos[1]< sizes[1];pos[1]++)
            for (pos[0]=0;pos[0]< sizes[0];pos[0]++)
                {
                if (grad != 0) {
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    for (int s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (1.0-mydata[i]/tmp)*tmpgrad[s];  // sum over all the pixels
                }
                else
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 

                if (res != 0)
                    res[i] = tmp; // save idiv image
                tmp=tmp-mydata[i]*log(tmp);  // fast version omitting constants
                if (resdiff != 0)
                    resdiff[i] = tmp; // save reduced idiff image: ATTENTION: This does not look like a good fit as the constant terms are omitted
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
                i++;
            }
    if (tmpgrad != 0)
        free(tmpgrad);
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
    double *data=0,tmp;  // Let's hope this one survives the individual calls
    double * params, * gparams,* fixedparams, result, * res=0, * resdiff=0, * grad=0;
    char * input_buf;
    int PsizeX,PsizeY, buflen; // GPsize, sizes[5],
    static double * mydata=0;
    static double * cparams=0;
    static int allocatedCparams=0;
    static int sizes[100],nd=1,totalsize=1,numdims;   // To hell with ppl who use more than 100 dimensions!
    static int docenter=0, switchxy=0;   // Boolean variable 
    enum method {mse, idiv,fidiv};
    static method mymethod=mse;
    
    if (nrhs != 1 && (nrhs < 3 || nrhs > 5)) 
    {
        printf ("Nr of parameters: %d\n",nrhs);
         mexErrMsgTxt("1 or 3-5 inputs required");
    }
    //if (mxIsChar(prhs[0]) != 1)
    //     mexErrMsgTxt("Input must be a string.");


    if (nrhs > 1)  // archieve the data by saving the pointer and remember which method to use
    {
    if (nrhs < 2)  
            mexErrMsgTxt("When submitting data, at least three arguments are required: data, method-string, dimensions!");
    if (nrhs > 5)  
            mexErrMsgTxt("When submitting data, at max 5 arguments are required: data, method-string, dimensions, docenter, switchXY!");
        
    /* Get the length of the input string. */
        const int *sz;
        sz = mxGetDimensions(prhs[0]);
        nd = mxGetNumberOfDimensions(prhs[0]);
        totalsize=1;
        for (int ii=0;ii<100;ii++)
        {
            sizes[ii]=1;
        }

        for (int ii=0;ii<nd;ii++)
        {
           sizes[ii]=sz[ii];  // save these values
           // printf("dim %d, size %d\n",ii,sizes[ii]);
           totalsize *= sz[ii];
        }
        if (nrhs > 4) {  // switchxy can be set
            switchxy=(int) (* mxGetPr(prhs[4]));
        }

        data = mxGetPr(prhs[0]);
        if (mydata != 0) free(mydata);
        mydata=(double *) calloc(totalsize,sizeof(double));
        if (switchxy != 0 && (nd > 1))
        {
            int sizexy=sizes[0]*sizes[1];
            int sizerest=totalsize/sizexy,tmp;
            for (int i=0;i<sizerest;i++) 
                for (int y=0;y<sizes[1];y++) 
                    for (int x=0;x<sizes[0];x++) 
                        mydata[i*sizexy+y+x*sizes[1]]=data[i*sizexy+x+y*sizes[0]];       // makes a copy to be able to access it in later calls
            tmp=sizes[0];
            sizes[0]=sizes[1];
            sizes[1]=tmp;
        }
        else
            for (int i=0;i<totalsize;i++) 
                mydata[i]=data[i];       // makes a copy to be able to access it in later calls

        if (nlhs != 0)
            mexErrMsgTxt("When submitting data, no output is returned!");
        // printf("Bufflen is %d x %d, pointer is %x, copied to %x\n",dataSizeX,dataSizeY,data,mydata);
    /* Get the length of the input string. */
    buflen = (mxGetM(prhs[1]) * mxGetN(prhs[1])) + 1;
    
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

    numdims=(int) (* mxGetPr(prhs[2]));  // get number of dimensions

    if (nrhs > 3) {  // docenter can be set
    docenter=(int) (* mxGetPr(prhs[3]));
    }


    mxFree(input_buf);
    }
    else // a single argument : Not the funtion has to be called and data simulated
    {
        if (mydata == 0)
            mexErrMsgTxt("No data was yet submitted. Cannot call function.");
        PsizeX = mxGetM(prhs[0]);
        if (mxGetN(prhs[0]) > 1)
            mexErrMsgTxt("All parameters should be in a single vector.");
        if ((PsizeX-1) % (numdims*2 +1)!=0)
            mexErrMsgTxt("Number of parameters (1 global and positions and sigmas) does not match with number of dimensions to fit.");
        else
        {
            PsizeY = (PsizeX-1)/(2*numdims+1);
            PsizeX = 2*numdims+1;  // Intensity, positions and variances
        }
        gparams = mxGetPr(prhs[0]);
        params = gparams+1;  // omit the one global (offset)
    	if (allocatedCparams< PsizeX*PsizeY+1 && docenter)
        {
            if (cparams != 0)
                free(cparams);
            cparams=(double *) calloc(PsizeX*PsizeY+1+100,sizeof(double));
            allocatedCparams=PsizeX*PsizeY+1+100;
        }

        //printf("ParamsizeX: %d, Y %d\n",PsizeX,PsizeY);
        //printf("dataSizeX: %d, Y %d\n",dataSizeX,dataSizeY);
        if (mydata != 0) {

            grad=0;
            if (nlhs >= 2)  // gradient is requested
            {
                int mysize[2];
                mysize[0]=PsizeX*PsizeY+1;
                mysize[1]=1;
                plhs[1] = mxCreateNumericArray(1, mysize, mxDOUBLE_CLASS, mxREAL);   
                if (mxGetM(plhs[1]) != mxGetM(prhs[0]))
                    mexErrMsgTxt("Gradient size has to be equal to parameter size");                    
                grad = mxGetPr(plhs[1]);
            }
            else 
                grad=0;
            
            if (nlhs >= 3)
            {
                plhs[2] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                res = mxGetPr(plhs[2]);
            }
            if (nlhs >= 4)
            {
                plhs[3] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                resdiff = mxGetPr(plhs[3]);
            }

            if (docenter)
                for (int d=0;d<PsizeY;d++)   // to account for center of image corresponding to zero
                {
                    for (int n=0;n<PsizeX;n++)   // to account for center of image corresponding to  zero
                        if (n>0 && n <= numdims)
                            cparams[d*PsizeX+n] = params[d*PsizeX+n] + floor((double) sizes[n-1]/2);  
                        else
                            cparams[d*PsizeX+n] = params[d*PsizeX+n];  
                }
            else
                cparams=params;
                
            //for (int y=0;y<PsizeY;y++) 
            //  for (int x=0;x<PsizeX;x++) 
            //    printf("Param %d, %d : %g\n",x,y,params[x+y*PsizeX]);

            switch (mymethod)
            {
                case mse:
                    result=do_mse(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                case idiv:
                    result=do_idiv(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                case fidiv:
                    result=do_fidiv(mydata, sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                default:
                    mexErrMsgTxt("Undefined method. Valid methods are : 'mse' and 'idiv'");
            }
                if (nlhs >= 1)
            {
                plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
                double * dp = mxGetPr(plhs[0]);
                (* dp) = result;
            }
            // printf("Result %g\n",result);
        }
        else
            mexErrMsgTxt("Please provide just the data matrix first");
    }
    
    return;
}
