#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>
#include <string.h>

#define SIZE       8
#define PI         3.141592653589793
#define TWO_PI     (2.0 * PI)
#define SWAP(a,b)  tempr=(a);(a)=(b);(b)=tempr

int nextPowerOf2(int m)
{
    int x = 1;
    while(x < m){
        x <<= 1;
    }
    return x;
}

double * scaleInputs(int * buf, int nItems, int nextPower)
{
    float * tempFloatArray;
    tempFloatArray = (float *) malloc(nextPower*sizeof(float));

    int fourthousands = 0;
    for( int i = 0, j = 0; i < nextPower; i++, j+= 2)
    {
        if(buf[i] > 0)
        {
            tempFloatArray[j] = (float)buf[i]/(float)2147483647.0;
            tempFloatArray[j+1] = (float)0.0;
            /*if(i == 4000)
            {
                fourthousands = j;
                printf("tempFloatArray[4000]: %f\n", tempFloatArray[j]);
                printf("tempFloatArray[4000+1]: %f\n", tempFloatArray[j+1]);
            }*/
        }
        else if(buf[i] < 0)
        {
            tempFloatArray[j] = (float)buf[i]/(float)2147483648.0;
            tempFloatArray[j+1] = (float)0.0;
            /*if(i == 4000)
            {
                fourthousands = j;
                printf("tempFloatArray[4000]: %f\n", tempFloatArray[j]);
                printf("tempFloatArray[4000+1]: %f\n", tempFloatArray[j+1]);
            }*/
        }
        else
        {
            tempFloatArray[j] = (float)0.0;
            tempFloatArray[j+1] = (float)0.0;
            /*if(i == 4000)
            {
                fourthousands = j;
                printf("tempFloatArray[4000]: %f\n", tempFloatArray[j]);
                printf("tempFloatArray[4000+1]: %f\n", tempFloatArray[j+1]);
            }*/
        }
    }

    //printf("j value of 4000: %d \n", fourthousands);
    //printf("tempFloatArray[8000]: %f\n", tempFloatArray[8000]);

    //printf("-------scale inputs--------\n");
    double * newBuf;
    int newSize = nextPower << 1;
    newBuf = (double *) malloc(newSize*sizeof(double));

    for(int i = 0; i < newSize; i++)
    {
        if( i < nextPower)
        {
            newBuf[i] = (double)(tempFloatArray[i]);
        }
        else
        {
            newBuf[i] = (double)0.0;
        }
    }

    //printf("newBuf[8000]: %lf\n", newBuf[8000]);
    return newBuf;
}


int * scaleOutputs(double * buf, int nItems)
{
    double posMax = 0;
    double negMax = 0;
    int * outputbuf;
    double * temp;

    outputbuf = (int *) malloc(nItems*sizeof(int));
    temp = (double *) malloc(nItems*sizeof(double));

    for(int i = 0; i < nItems; i++)
    {
        if(buf[i] > 0 && buf[i] > posMax)
        {
            posMax = buf[i];
        }
        else if(buf[i] < 0 && buf[i] < negMax)
        {
            negMax = buf[i];   
        }
    }

    //printf("posMax: %f\n", posMax);
    //printf("negMax: %f\n", negMax);
    //printf("buf[100]: %f\n", buf[100]);

    for(int i = 0; i < nItems; i++)
    {
        if(buf[i] > 0)
        {
            temp[i] = buf[i]/posMax;
        }
        else if(buf[i] < 0)
        {
            temp[i] = -1 * buf[i]/negMax;  
        }
    }

    //printf("temp[100]: %f\n", temp[100]);

    for(int i = 0; i < nItems; i++)
    {
        if(temp[i] > 0)
        {
            outputbuf[i] = (int)(temp[i] * (double)2147483647.0);
        }
        else if(temp[i] < 0)
        {
            outputbuf[i] = (int)(temp[i] * (double)2147483648.0);
        }
        else
        {
            outputbuf[i] = 0;
        }
    }

    //printf("outputbuf[100]: %d\n", outputbuf[100]);
    return outputbuf;
}

/*Reference Smith, p315 of FFT Convolution
        re Y[k] = re X[k] re H[k] - im X[k] im H[k]
        im Y[k] = im X[k] re H[k] + re X[k] im H[k] */

double *complexMultiply(double x[], double h[], int ii){
    double *y; 
    y = (double *) malloc(ii*sizeof(double));
    for(int i = 0; i < ii-2; i+= 2){
        double reXH = x[i] * h[i];      //re X[k] re H[k]
        double imXH = x[i+1] * h[i+1];  //im X[k] im H[k]
        double reY = reXH - imXH;       //re Y[k] = re X[k] re H[k] - im X[k] im H[k]
        y[i] = reY;                     //real Y

        double irXH = x[i+1] * h[i];    //im X[k] re H[k]
        double riXH = x[i] * h[i+1];    //re X[k] im H[k]
        double imY = irXH + riXH;       //im Y[k] = im X[k] re H[k] + re X[k] im H[k]
        y[i+1] = imY;                   //imaginary Y
    }
    return y;
}

//  The four1 FFT from Numerical Recipes in C,
//  p. 507 - 508.
//  Note:  changed float data types to double.
//  nn must be a power of 2, and use +1 for
//  isign for an FFT, and -1 for the Inverse FFT.
//  The data is complex, so the array size must be
//  nn*2. This code assumes the array starts
//  at index 1, not 0, so subtract 1 when
//  calling the routine (see main() below).

void four1(double data[], int nn, int isign)
{
    unsigned long n, mmax, m, j, istep, i;
    double wtemp, wr, wpr, wpi, wi, theta;
    double tempr, tempi;

    n = nn << 1;
    j = 1;

    for (i = 1; i < n; i += 2) {
    if (j > i) {
        SWAP(data[j], data[i]);
        SWAP(data[j+1], data[i+1]);
    }
    m = nn;
    while (m >= 2 && j > m) {
        j -= m;
        m >>= 1;
    }
    j += m;
    }

    mmax = 2;
    while (n > mmax) {
    istep = mmax << 1;
    theta = isign * (6.28318530717959 / mmax);
    wtemp = sin(0.5 * theta);
    wpr = -2.0 * wtemp * wtemp;
    wpi = sin(theta);
    wr = 1.0;
    wi = 0.0;
    for (m = 1; m < mmax; m += 2) {
        for (i = m; i <= n; i += istep) {
        j = i + mmax;
        tempr = wr * data[j] - wi * data[j+1];
        tempi = wr * data[j+1] + wi * data[j];
        data[j] = data[i] - tempr;
        data[j+1] = data[i+1] - tempi;
        data[i] += tempr;
        data[i+1] += tempi;
        }
        wr = (wtemp = wr) * wpr - wi * wpi + wr;
        wi = wi * wpr + wtemp * wpi + wi;
    }
    mmax = istep;
    }
}


int main(int argc, char * argv[])
{
    if(argc != 4)
    {
        printf("Please enter the 3 files in the form of /convolve input.wav impulse.wav output.wav\n");
        exit(-1);
    }

    char * inputFileName = argv[1];
    char * impulseFileName = argv[2];
    char * outputFileName = argv[3];
    //printf("here\n");

    SNDFILE *sf;
    SF_INFO info;
    int num_channels;
    int num, num_items;
    int *buf;
    int f,sr,c;
    
    /* Open the WAV file. */
    info.format = 0;
    sf = sf_open(inputFileName,SFM_READ,&info);
    if (sf == NULL)
    {
        printf("Failed to open the file.\n");
        exit(-1);
    }

    /* Print some of the info, and figure out how much data to read. */
    f = info.frames;
    sr = info.samplerate;
    c = info.channels;
    num_items = f*c;        // number of items need to read from input file
    printf("frames=%d\n",f);
    printf("samplerate=%d\n",sr);
    printf("channels=%d\n",c);
    //printf("format=%d\n",info.format);
    //printf("sections=%d\n",info.sections);
    //printf("seekable=%d\n",info.seekable);
    printf("num_items=%d\n",num_items);
    

    /* Allocate space for the data to be read, then read it. */
    buf = (int *) malloc(num_items*sizeof(int));
    num = sf_read_int(sf,buf,num_items);
    sf_close(sf);
    printf("Read %d items\n",num);

    printf("buf[4000]:%d\n", buf[4000] );


    printf("\n");


    printf("-----IMPULSE INFO -------------------------\n");

    SNDFILE *impulse;
    SF_INFO impulseInfo;
    int impulse_num_channels;
    int impulse_num, impulse_num_items;
    int *impulse_buf;
    int i_f,i_sr,i_c;
    
    /* Open the WAV file. */
    impulseInfo.format = 0;
    impulse = sf_open(impulseFileName,SFM_READ, &impulseInfo);
    if (impulse == NULL)
    {
        printf("Failed to open the impulse file.\n");
        exit(-1);
    }
    //printf("here2\n");

    /* Print some of the info, and figure out how much data to read. */
    i_f = impulseInfo.frames;
    i_sr = impulseInfo.samplerate;
    i_c = impulseInfo.channels;
    impulse_num_items = i_f*i_c;        // number of items need to read from input file
    printf("frames=%d\n",i_f);
    printf("samplerate=%d\n",i_sr);
    printf("channels=%d\n",i_c);
    //printf("format=%d\n",impulseInfo.format);
    //printf("sections=%d\n",impulseInfo.sections);
    //printf("seekable=%d\n",impulseInfo.seekable);
    printf("num_items=%d\n",impulse_num_items);
    

    /* Allocate space for the data to be read, then read it. */
    impulse_buf = (int *) malloc(impulse_num_items*sizeof(int));
    impulse_num = sf_read_int(impulse,impulse_buf,impulse_num_items);
    sf_close(impulse);
    printf("Read %d items\n",impulse_num);

    printf("impulse_buf[4000]:%d\n", impulse_buf[4000] );

    // find the maximum size
    int max;
    if(num >= impulse_num)
    {
        max = num;
    }
    else
    {
        max = impulse_num;
    }

    printf("----------------\n");
    printf("Max impulse : %d\n", max);
    int nextPower = nextPowerOf2(max);

    printf("----------------\n");
    printf("Next power : %d\n", nextPower);

    int * inputP2;
    inputP2 = (int *) malloc((nextPowerOf2(max))*sizeof(int));
    int * impulseP2;
    impulseP2 = (int *) malloc((nextPowerOf2(max))*sizeof(int));


    // jamming tuning technique
    for(int i = 0; i < (nextPowerOf2(max)); i++)
    {
        inputP2[i] = 0;
    }

    for(int i = 0; i < (nextPowerOf2(max)); i++)
    {
        impulseP2[i] = 0;
    }

    printf("before memory copy\n");
    printf("buf[4000]: %d\n", buf[4000]);
    printf("impulse_buf[4000]: %d\n", impulse_buf[4000] );

    // copy all memory from input buffers to the arrays
    memcpy(inputP2, &buf[0], num);
    memcpy(impulseP2, &impulse_buf[0], impulse_num);

    // check if the memory copy was successful
    printf("----------------\n");
    printf("inputP2[4000]:%d\n", inputP2[4000] );
    printf("impulseP2[4000]: %d\n", impulseP2[4000] );

    //scaling items in input buffer and impulse buffer
    printf("\n");
    printf("#### newInputBuf ####\n");
    double * newInputBuf;
    newInputBuf = (double *)malloc((nextPower << 1) * sizeof(double));
    newInputBuf = scaleInputs(inputP2, num_items, nextPower);
    printf("\n");
    printf("#### newImpulseBuf ####\n");
    double * newImpulseBuf;
    newImpulseBuf = (double *)malloc((nextPower << 1) * sizeof(double));
    newImpulseBuf = scaleInputs(impulseP2, impulse_num_items, nextPower);

    // check if the new double value was scaled correctly
    printf("----------------\n");
    printf("newInputBuf[8000]:%lf\n", newInputBuf[8000] );
    printf("newImpulseBuf[8000]: %lf\n", newImpulseBuf[8000] );


    four1(newInputBuf - 1, (nextPowerOf2(max)), 1);
    four1(newImpulseBuf - 1, (nextPowerOf2(max)), 1);

    double * multiply;
    multiply = (double *) malloc((nextPower<<1)*sizeof(double));
    multiply = complexMultiply(newInputBuf, newImpulseBuf, nextPower);

    four1(multiply - 1, nextPower, -1);

    //scale after doing ifft
    for(int i = 0; i < nextPower << 1; i++){
        multiply[i] /= (double) nextPower;
    }


    /*printf("-----------------------\n");

    printf("convolving\n");
    //convolve the current one
    convolve(newInputBuf, num, newImpulseBuf, impulse_num, output_buf, output_buf_size);
    printf("finished convolving\n");


    //int temp_output_size = 19999;
    int * final_output_buf = scaleOutputs(output_buf, output_buf_size);
    printf("here1\n");*/
    //float *output_buf;
    int *outputvalues;
    outputvalues = (int *) malloc(nextPower*sizeof(int));

    // tuning: eliminate common subexpressions
    int output_buf_size = num + impulse_num - 1;

    //allocate space for float data in output data;
    //output_buf = (float *) malloc(output_buf_size*sizeof(float));

    outputvalues = scaleOutputs(multiply, num + impulse_num - 1);


    SNDFILE *out;
    SF_INFO outInfo;

    outInfo.format = 0;
    outInfo.frames = f;
    outInfo.samplerate = 44100;
    outInfo.channels = 1;
    outInfo.format = SF_FORMAT_WAV | SF_FORMAT_PCM_16 | SF_ENDIAN_LITTLE;

    out = sf_open(outputFileName, SFM_WRITE, &outInfo);

    if (out == NULL)
    {
        printf("cannot open file to write\n");
        exit(-1);
    }

    printf("here2\n");

    int outputSize = sf_write_int(out, outputvalues, num + impulse_num - 1); // temp_output_size should be actual output size
    
    printf("here3\n");
    sf_close(out);
    return 0;
}