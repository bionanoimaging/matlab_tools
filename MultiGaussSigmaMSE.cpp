/* Compares a number of simulated Gaussian Spots with experimental data

 * can account for a ROI around every spot to make it fast
 * Examle : a=exp(-((xx(20,20)-2).^2+(yy(20,20)-1.2).^2)/20)
           [fitted,params]=FitDataNDFast([0 1 0; 1 0 0],a)
 */

#include "mex.h"
#include <stdio.h>
#include <string.h>
#include <math.h>

static int sizexy;  
static int sizey;  
static int sizerest;
static int totalsize=1;

static int pixelPos=0; // running variables: current pixel number
static int prevPixelPos=0; // running variables: previous pixel number
static int particleNum=0; // current particle
static int listPos=0;     // current position in list
static int listLength=0;
static int * pixellist=0;  // holds the information (only for the masked pixels) about: pixelnumber and spotindex, from which the position can be calculated
static double * listval=0;   // list of corresponding values
static int NumParticles=0;  // number of particles 
static double * mydata=0;
static double * tmpgrad=0;  // temporary gradient


void CurrentPos(int * pos) // calculates the position given a pixel number
{
    int tmp;
    pos[2]=pixelPos/sizexy;
    tmp=pixelPos%sizexy;
    pos[1]=tmp/sizey;
    pos[0]=tmp%sizey;    
    // printf("pixelPos : %d, pos: %d %d %d particleNum %d\n",pixelPos,pos[0],pos[1],pos[2],particleNum);
}

double pixelVal()  // returns 1 if the current is the same as the previous particle, 0 if it was the last particle for this pixel and -1 if there are no more pixels
{
    if (pixellist)
    {
        // return listPos;
        return listval[listPos]; // value of current pixel
    }
    else  // use ordinary indexing
    {
        return mydata[pixelPos];
    }    
}

bool isValidPos()
{
    if (pixellist)
        if (listPos >= listLength)
            return false;
        else
            return true;
    else
        if (pixelPos >= totalsize)
            return false;
        else
            return true;        
}

int nextPos()  // returns 1 if the current is the same as the previous particle, 0 if it was the last particle for this pixel and -1 if there are no more pixels
{
    if (pixellist)
    {
        int oldParticleNum=particleNum;
        listPos++;  //global pixel index
        if (listPos > listLength)
            return -1;
        particleNum=pixellist[listPos*2];
        prevPixelPos=pixelPos;
        pixelPos=pixellist[listPos*2+1];
        // printf("OK pixelPos : %d, particleNum %d\n",pixelPos,particleNum);
        if (pixelPos == prevPixelPos) // is it still the same pixel, but simply a new particle
            return 1;
        else
            return 0;  // this is a new particle
    }
    else  // use ordinary indexing
    {
        if (pixelPos > totalsize)
            return -1;
        if (particleNum<NumParticles-1)
        {
            // printf("OK pixelPos : %d, particleNum %d\n",pixelPos,particleNum);
            particleNum++;
            return 1;
        }
        else
        {
            // printf("pixelPos : %d, particleNum %d\n",pixelPos,particleNum);
            particleNum=0;
            prevPixelPos=pixelPos;
            pixelPos++;
            // if (pixelPos>10)
            //    mexErrMsgTxt("Check 3b");
            return 0;
        }
    }    
}

// The function below calculates the N-dimensional Gaussian value at a given position pos
double dosim(double * params, int paramSize, int numparams, double * gparams, int * pos)
{
    double result=gparams[0],sq,ssq;  // offset is accounted for first gparams[0]
    int idx;
    //printf("PosX : %d, y %d\n",pos[0],pos[1]);
    do {
        int numdim=(paramSize-1)/2,d;
        ssq=0;
        for (d=0;d<numdim;d++)
        {
            idx=paramSize*particleNum+d+1;
            sq=(pos[d]-params[idx])/params[idx+numdim];  // positions are in [paramSize*n+d+1], sigmas in [paramSize*n+numdim+d+1]
        //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
         ssq +=sq*sq;}
        result += params[paramSize*particleNum]*exp(-ssq);                              // intensities are in params[paramSize*n]
        //printf("Intensity : %g\n",params[paramSize*n]);
    }  while (nextPos()==1); // this is the same pixel position, but differen particles.

    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    return result;
}

// The function below calculates the Value (return value) and the Gradient of a pixel at position pos
double doGradSim(double * mygradient, double * params, int paramSize, int numparams, double * gparams, int * pos)
{
    double result=0.0,spotresult,sq,ssq;  
    // printf("PosX : %d, y %d\n",pos[0],pos[1]);
    do {  // iterates over all particles of the same pixel, gradient results are stored in a temporary gradient
        int numdim=(paramSize-1)/2;   // dimension dependent are intensities and sigmas
        int idx,d;
        ssq=0;
        for (d=0;d<numdim;d++) {
            idx=paramSize*particleNum+d+1;
            // printf("dimension d %d / %d, idx %d\n",d,numdim,idx);
            sq=(pos[d]-params[idx])/params[idx+numdim];  // positions are in [paramSize*n+d+1], sigmas in [paramSize*n+numdim+d+1]
            //printf("Param %d %d: %g\n",d,n,params[paramSize*n+d+1]);
            ssq +=sq*sq;
            }
        spotresult = exp(-ssq);                                        // intensities are in params[paramSize*n]
        mygradient[1+paramSize*particleNum] = spotresult;              // gradient with respect to I
        spotresult *= params[paramSize*particleNum];                   // mulitply by I
        result += spotresult;
        for (d=0;d<numdim;d++) {
            idx=paramSize*particleNum+d+1;
            sq= 2.0*(pos[d]-params[idx])/(params[idx+numdim] * params[idx+numdim]);
            mygradient[1+idx] = spotresult * sq ;   // positional gradient
            mygradient[1+idx+numdim] = spotresult * sq * (pos[d]-params[idx])/params[idx+numdim];   // sigma gradient
            // printf("dimension d %d / %d, idx %d, %d\n",d,numdim,idx,idx+numdim);
        }
        //printf("Intensity : %g\n",params[paramSize*n]);
    }    while (nextPos()==1); // this is the same pixel
    result += gparams[0];               // offset is accounted for first gparams[0]
    // gradient with respect to global background is always 1.0
    mygradient[0] = 1.0;
    //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
    // mexErrMsgTxt("Check GradDone");
    return result;
}

// this function computes the meas square error comparing data with simulation
// params is the array of spot positions
// paramsg are the global parameters
double do_mse(int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double *resdiff, double *grad)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    //for (int y=0;y<numparams;y++) 
    //   for (int x=0;x<paramsize;x++) 
    //      printf("MSE: Param %d, %d : %g\n",x,y,params[x+y*paramsize]);
    double result=0,tmp,tmp2,* tmpgrad=0,val;
    int mypos[3],* pos=mypos,s;
    if (grad != 0)
    {
        if (tmpgrad) {free(tmpgrad);tmpgrad=0;}
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    pos[0]=0;pos[1]=0;pos[2]=0;
    while (isValidPos())  // the inner loop in doGradSim or dosim will update the pixel and list postions
                {
                val=pixelVal();
                CurrentPos(pos);  // update position information
                if (grad != 0) {
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    tmp2 = val-tmp;
                    for (s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (-2.0)*tmpgrad[s]*tmp2;  // sum over all the pixels
                }
                else
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 
                if (res != 0)
                    res[prevPixelPos] = tmp; // save simulation
                tmp=val-tmp;
                if (resdiff != 0)
                    resdiff[prevPixelPos] = tmp; // save difference                
                result += tmp*tmp;     // sum of squares over all pixels
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
            }
    if (tmpgrad != 0)
        free(tmpgrad);
    return result;
}

double do_idiv(int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff, double *grad)  // int numgparams, 
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp,val,* tmpgrad=0;
    int pos[3],s;
    // mexErrMsgTxt("Check 1");
    if (grad != 0)
    {
        if (tmpgrad) {free(tmpgrad);tmpgrad=0;}
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    while (isValidPos())  // the inner loop in doGradSim or dosim will update the pixel and list postions
       {
          val=pixelVal();
          CurrentPos(pos);  // update position information
                if (grad != 0) {
                    // mexErrMsgTxt("Check 18a");
                    // printf("tmp : %g, val : %g, Pos: %d, Length %d\n",tmp,val,listPos,listLength);
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    for (s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (1.0-val/tmp)*tmpgrad[s];  // sum over all the pixels
                    // printf("tmp : %g, GradX: %g, SumGrad %g, PosX: %d, val : %g, Pos: %d, Length %d\n",tmp,tmpgrad[1],grad[1],pos[0],val,listPos,listLength);
                }
                else
                {
                    // mexErrMsgTxt("Check 18b");
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 
                }
          
                if (res != 0)
                    res[prevPixelPos] = tmp; // save idiv image
                if (val !=0)
                    tmp=val*log(val/tmp)-(val-tmp);  // i-divergence with sterling's approximation
                else
                    tmp=-(val-tmp);
                if (resdiff != 0)
                    resdiff[prevPixelPos] = tmp; // save difference                
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
            }
    if (tmpgrad != 0)
        free(tmpgrad);
    //printf("result : %g, Params = %g",result,gparams[0]); for (int s=1;s<paramsize*numparams+1;s++) printf(", %g",gparams[s]); printf("\n");
    //printf("result : %g, Grad = %g",result,grad[0]); for (int s=1;s<paramsize*numparams+1;s++) printf(", %g",grad[s]); printf("\n");
    // mexErrMsgTxt("Check 4");
    return result;
}

double do_fidiv(int * sizes, double * params, int paramsize, int numparams, double * gparams, double * res, double * resdiff, double *grad)
{
    //printf("Bufflen is %d x %d, pointer is %x\n",sizes[0],sizes[1],mydata);
    double result=0,tmp,*tmpgrad=0,val;
    int i=0,pos[3],s;
    if (grad != 0)
    {
        if (tmpgrad) {free(tmpgrad);tmpgrad=0;}
        tmpgrad=(double *) calloc(paramsize*numparams+1,sizeof(double));
        for (s=0;s<paramsize*numparams+1;s++)
            grad[s] = 0.0;
    }    
    while (isValidPos())  // the inner loop in doGradSim or dosim will update the pixel and list postions
       {
          val=pixelVal();
          CurrentPos(pos);  // update position information
                if (grad != 0) {
                    tmp = doGradSim(tmpgrad, params, paramsize, numparams, gparams, pos);
                    for (s=0;s<paramsize*numparams+1;s++)
                            grad[s] += (1.0-val/tmp)*tmpgrad[s];  // sum over all the pixels
                }
                else
                    tmp = dosim(params, paramsize, numparams, gparams, pos); 

                if (res != 0)
                    res[prevPixelPos] = tmp; // save idiv image
                tmp=tmp-val*log(tmp);  // fast version omitting constants
                if (resdiff != 0)
                    resdiff[prevPixelPos] = tmp; // save reduced idiff image: ATTENTION: This does not look like a good fit as the constant terms are omitted
                result += tmp;
                //printf("PosX : %d, y %d, val= %g\n",pos[0],pos[1],result);
                //printf("tmp = %g\n",tmp);
                //printf("mydata = %g\n",mydata[i]);
            }
    if (tmpgrad != 0)
        free(tmpgrad);
    return result;
}


bool inROI1D(int * pos,double * params, double * ROI)
{
    return (pos[0]>params[0]-ROI[0]) && (pos[0]<params[0]+ROI[0]);
}

bool inROI2D(int * pos,double * params, double * ROI)
{
    return (pos[0]>params[0]-ROI[0]) && (pos[0]<params[0]+ROI[0]) &&
            (pos[1]>params[1]-ROI[1]) && (pos[1]<params[1]+ROI[1]);
}

bool inROI3D(int * pos,double * params, double * ROI)
{
    return (pos[0]>params[0]-ROI[0]) && (pos[0]<params[0]+ROI[0]) &&
            (pos[1]>params[1]-ROI[1]) && (pos[1]<params[1]+ROI[1]) &&
            (pos[2]>params[2]-ROI[2]) && (pos[2]<params[2]+ROI[2]);
}

bool inROI(int numdims, int * pos,double * params, double * ROI)
{
    switch (numdims)
    {
            case 1:
                return inROI1D(pos,params,ROI);
            case 2:
                return inROI2D(pos,params,ROI);
            case 3:
                return inROI3D(pos,params,ROI);
            default:
                mexErrMsgTxt("ROISize unsupported number of dimensions");
    }
    return false;
}


                                
void mexFunction( int nlhs, mxArray *plhs[],
                  int nrhs, const mxArray *prhs[] )
{
   //printf ("Globals: sizey %d size xy %d, sizerest %d, mydata %d\n",sizey, sizexy, sizerest, mydata);
   //printf ("Globals: pixelPos %d ParticleNum %d, listPos %d, listLength %d\n",pixelPos, particleNum, listPos, listLength);
   //printf ("Globals: pixellist %d listval %d, NumParticles %d\n",pixellist, listval, NumParticles);

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
    int pos[3],d,n;
    static double * cparams=0;
    static int allocatedCparams=0;
    static int sizes[100],nd=1,numdims;   // To hell with ppl who use more than 100 dimensions!
    static int docenter=0, switchxy=0;   // Boolean variable 
    enum method {mse, idiv,fidiv};
    static enum method mymethod=mse;
    
    if (nrhs != 1 && (nrhs < 3 || nrhs > 7)) 
    {
        printf ("Nr of parameters: %d\n",nrhs);
         mexErrMsgTxt("1, 3-7 inputs required");  // 1: SpotList only, 3-5: Data and other parameters, 6-7 spot positions and ROI sizes
    }
    //if (mxIsChar(prhs[0]) != 1)
    //     mexErrMsgTxt("Input must be a string.");


    if (nrhs > 1)  // archieve the data by saving the pointer and remember which method to use
    {
        const int *sz;
        int ii,n;
        
    if (nrhs < 2)  
            mexErrMsgTxt("When submitting data, at least three arguments are required: data, method-string, dimensions!");
    if (nrhs > 8)
            mexErrMsgTxt("When submitting data, at max 5 arguments are required: data, method-string, dimensions, docenter, switchXY, initialSpotList, ROIsizes!");
        
    /* Get the length of the input string. */
        sz = mxGetDimensions(prhs[0]);
        nd = mxGetNumberOfDimensions(prhs[0]);
        totalsize=1;
        for (ii=0;ii<100;ii++)
        {
            sizes[ii]=1;
        }

        for (ii=0;ii<nd;ii++)
        {
           sizes[ii]=sz[ii];  // save these values
           // printf("dim %d, size %d\n",ii,sizes[ii]);
           totalsize *= sz[ii];
        }
        if (nrhs > 4) {  // switchxy can be set
            switchxy=(int) (* mxGetPr(prhs[4]));
        }

        data = mxGetPr(prhs[0]);

        listLength=0;
        if (pixellist != 0) free(pixellist);
        if (listval != 0) free(listval);
        pixellist=0;listval=0;

        numdims=(int) (* mxGetPr(prhs[2]));  // get number of dimensions

        sizey=sizes[1];
        sizexy=sizes[0]*sizes[1];
        sizerest=totalsize/sizexy;
        if (nrhs > 5) {  // SpotList to use to build a pixelmask
            double ROIs[3]={5.0, 5.0, 5.0}, * ROI;
            gparams = mxGetPr(prhs[5]);
            printf ("Size of argument 5 is : % d, numdims is %d\n",mxGetM(prhs[5]),numdims);
            PsizeX = 2*numdims+1;  // Intensity, positions and variances
            PsizeY = (mxGetM(prhs[5])-1)/PsizeX;
            if (PsizeY * (2*numdims+1) != mxGetM(prhs[5])-1)
                mexErrMsgTxt("Error: Wrong number of initial arguments for this number of dimensions. Arguments should include global background");
                
            NumParticles=PsizeY;  // will never be changed from now on
            ROI=ROIs;
            if (nrhs > 6)  // ROI sizes
                {
                    ROI = mxGetPr(prhs[6]);
                    if (mxGetN(prhs[6]) != numdims)
                        printf ("WARNING: Wrong size of ROIs, has to correspond to number of dimensions. Using only valid ones.\n");
                        // mexErrMsgTxt("Error: Wrong size of ROIs, has to correspond to number of dimensions");
                }
            listLength=PsizeY;
            for (n=0;n<numdims;n++) listLength *= (int) (2*ROI[n]+1);  // stores the length of the list
            listLength++;
            // printf ("Allocating pixellist length %d for %d particles, PsizeX %d\n",listLength,PsizeY,PsizeX);
            if (pixellist) {free(pixellist);pixellist=0;}
            if (listval) {free(listval); listval=0;}
            pixellist=(int *) calloc(listLength*2, sizeof(int));  // possibly too big, depending on clipped pixels
            listval=(double *) calloc(listLength, sizeof(double));  // possibly too big, depending on clipped pixels
            params = gparams+1;  // omit the one global (offset)
            listPos=0;
            for (pos[2]=0;pos[2]<sizes[2];pos[2]++)  // iterate over all pixels
                for (pos[1]=0;pos[1]<sizes[1];pos[1]++)
                    for (pos[0]=0;pos[0]<sizes[0];pos[0]++) {
                        for (n=0;n<PsizeY;n++)  // go through all the particles
                            if (inROI(numdims,pos, &params[1+PsizeX*n], ROI)) { // only pixels in the ROI of a particular spot are stored in the list
                                pixellist[2*listPos]=n;   // store the particle number
                                pixelPos=pos[2]*sizexy+pos[1]*sizes[0]+pos[0];
                                pixellist[2*listPos+1]=pixelPos;  // store the pixel number
                                if (switchxy != 0 && (nd > 1))
                                    pixelPos=pos[2]*sizexy+pos[0]*sizey+pos[1];  // store the pixel number
                                listval[listPos]=data[pixelPos]; // store the data value of this pixel. Replicated for each spot
                                // printf ("Filling pixel %d, particle %d, listPos %d, value %g \n",pixelPos,n,listPos,data[pixelPos]);
                                listPos++;
                                }
                    }
            listLength=listPos;  // stores the length of the list
            listPos=0;
//            mexErrMsgTxt("Check 1");
            }
        else // no spotlist. The data has to be copied
            {
            int i;
                // printf("reallocating mydata\n");
                if (mydata != 0) free(mydata);
                mydata=(double *) calloc(totalsize, sizeof(double));
                    
                if (switchxy != 0 && (nd > 1))
                 {
                        int tmp,i,x,y;
                        for (i=0;i<sizerest;i++)
                            for (y=0;y<sizes[1];y++)
                                for (x=0;x<sizes[0];x++)
                                    mydata[i*sizexy+y+x*sizes[1]]=data[i*sizexy+x+y*sizes[0]];       // makes a copy to be able to access it in later calls
                        tmp=sizes[0];
                        sizes[0]=sizes[1];
                        sizes[1]=tmp;
                }
                else
                    for (i=0;i<totalsize;i++)
                        mydata[i]=data[i];       // makes a copy to be able to access it in later calls
        }

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

    if (nrhs > 3) {  // docenter can be set
    docenter=(int) (* mxGetPr(prhs[3]));
    }


    mxFree(input_buf);
    }
    else // a single argument : Now the funtion has to be called, data simulated and the error computed according to the chosen norm
    {
        if (mydata == 0 && (pixellist == 0 || listLength ==0))
            mexErrMsgTxt("No data was yet submitted. Cannot call function.");
        PsizeX = mxGetM(prhs[0]);
        if (mxGetN(prhs[0]) > 1)
            mexErrMsgTxt("All parameters should be in a single vector.");
        if ((PsizeX-1) % (numdims*2 +1)!=0)
            mexErrMsgTxt("Number of parameters (1 global and positions and sigmas) does not match with number of dimensions to fit.");
        else
        {
            PsizeY = (PsizeX-1)/(2*numdims+1);
            if (pixellist && (NumParticles != PsizeY))
            {
                printf("NumParticles %d, Stored NumParticles %d\n",PsizeY,NumParticles);
                mexErrMsgTxt("Number of particles to fit cannot be changed after initialization when using ROIs.");
            }   
            NumParticles = PsizeY; // will be used in the inner loop
            PsizeX = 2*numdims+1;  // Intensity, positions and variances
        }
        gparams = mxGetPr(prhs[0]);
        params = gparams+1;  // omit the one global (offset)
    	if (allocatedCparams< PsizeX*PsizeY+1 && docenter)
        {
            if (cparams != 0) {free(cparams);cparams=0;}
            cparams=(double *) calloc(PsizeX*PsizeY+1+100,sizeof(double));
            allocatedCparams=PsizeX*PsizeY+1+100;
        }

        //printf("ParamsizeX: %d, Y %d\n",PsizeX,PsizeY);
        //printf("dataSizeX: %d, Y %d\n",dataSizeX,dataSizeY);
        if (mydata != 0 || pixellist != 0) {

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
            
            if (nlhs >= 3)  // generates a result
            {
                int p;
                plhs[2] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                res = mxGetPr(plhs[2]);
                if (pixellist)
                    for (p=0;p<totalsize;p++) res[p]=gparams[0];   // fill in the Background value
            }
            if (nlhs >= 4)  // generate the difference image
            {
                plhs[3] = mxCreateNumericArray(nd, sizes, mxDOUBLE_CLASS, mxREAL);
                //plhs[1] = mxCreateDoubleMatrix(sizes,nd, mxREAL);
                resdiff = mxGetPr(plhs[3]);
            }

            if (docenter)
                for (d=0;d<PsizeY;d++)   // to account for center of image corresponding to zero
                {
                    for (n=0;n<PsizeX;n++)   // to account for center of image corresponding to  zero
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

            particleNum=0;
            listPos=0;
            prevPixelPos=-1;
            if (pixellist) {
                particleNum=pixellist[2*listPos]; // running variables
                pixelPos=pixellist[2*listPos+1]; // running variables
            }
            else
            {
                particleNum=0;
                pixelPos=0;
            }
            switch (mymethod)
            {
                case mse:
                    result=do_mse(sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                case idiv:
                    result=do_idiv(sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                case fidiv:
                    result=do_fidiv(sizes, cparams, PsizeX, PsizeY, gparams, res, resdiff, grad); // GPsize,
                break;
                default:
                    mexErrMsgTxt("Undefined method. Valid methods are : 'mse' and 'idiv'");
            }
                if (nlhs >= 1)
            {
                double * dp;
                plhs[0] = mxCreateDoubleMatrix(1,1, mxREAL);
                dp = mxGetPr(plhs[0]);
                (* dp) = result;
            }
            // printf("Result %g\n",result);
        }
        else
            mexErrMsgTxt("Please provide just the data matrix first");
    }
    
    return;
}
