#include <stdio.h>
#include <stdlib.h>
#include <sndfile.h>
#include <math.h>

float * scaleInputs(int * buf, int nItems)
{
    printf("-------scale inputs--------\n");
    float * newBuf;
    newBuf = (float *) malloc(nItems*sizeof(float));
    printf(" buf[5]: %d\n", buf[5]);

    for(int i = 0; i < nItems; i++)
    {
        if(buf[i] > 0)
        {
            newBuf[i] = (float)buf[i]/2147483647.0f;
        }
        else if(buf[i] < 0)
        {
            newBuf[i] = (float)buf[i]/2147483648.0f;
        }
        else
        {
            newBuf[i] = 0.0f;
        }
    }

    printf("newBuf[5]: %f\n", newBuf[5]);
    return newBuf;
}

void convolve(float x[], int N, float h[], int M, float y[], int P)
{
    /*  Make sure the output buffer is the right size: P = N + M - 1  */
    if (P != (N + M - 1)) 
    {
        printf("Output signal vector is the wrong size\n");
        printf("It is %-d, but should be %-d\n", P, (N + M - 1));
        printf("Aborting convolution\n");
        return;
    }

    /*  Clear the output buffer y[] to all zero values  */
    for (int n = 0; n < P; n++)
    {   
        //printf("here\n");
        y[n] = 0.0;
    }
    

    /*  Do the convolution  */
    /*  Outer loop:  process each input value x[n] in turn  */
    for (int n = 0; n < N; n++) 
    {
        /*  Inner loop:  process x[n] with each sample of h[]  */
        for (int m = 0; m < M; m++)
            y[n+m] += x[n] * h[m];
    }
}

int * scaleOutputs(float * buf, int nItems)
{
    float posMax = 0;
    float negMax = 0;
    int * outputbuf;
    float * temp;

    outputbuf = (int *) malloc(nItems*sizeof(int));
    temp = (float *) malloc(nItems*sizeof(float));

    // find the max positive number and min negative number for scaling 
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

    printf("posMax: %f\n", posMax);
    printf("negMax: %f\n", negMax);
    printf("buf[100]: %f\n", buf[100]);

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

    printf("temp[100]: %f\n", temp[100]);

    for(int i = 0; i < nItems; i++)
    {
        if(temp[i] > 0)
        {
            outputbuf[i] = (int)(temp[i] * 2147483647.0f);
        }
        else if(temp[i] < 0)
        {
            outputbuf[i] = (int)(temp[i] * 2147483648.0f);
        }
        else
        {
            outputbuf[i] = 0;
        }
    }

    printf("outputbuf[100]: %d\n", outputbuf[100]);
    return outputbuf;
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
    printf("format=%d\n",info.format);
    printf("sections=%d\n",info.sections);
    printf("seekable=%d\n",info.seekable);
    printf("num_items=%d\n",num_items);
    

    /* Allocate space for the data to be read, then read it. */
    buf = (int *) malloc(num_items*sizeof(int));
    num = sf_read_int(sf,buf,num_items);
    sf_close(sf);
    printf("Read %d items\n\n",num);

    

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

    /* Print some of the info, and figure out how much data to read. */
    i_f = impulseInfo.frames;
    i_sr = impulseInfo.samplerate;
    i_c = impulseInfo.channels;
    impulse_num_items = i_f*i_c;        // number of items need to read from input file
    printf("frames=%d\n",i_f);
    printf("samplerate=%d\n",i_sr);
    printf("channels=%d\n",i_c);
    printf("format=%d\n",impulseInfo.format);
    printf("sections=%d\n",impulseInfo.sections);
    printf("seekable=%d\n",impulseInfo.seekable);
    printf("num_items=%d\n",impulse_num_items);
    

    /* Allocate space for the data to be read, then read it. */
    impulse_buf = (int *) malloc(impulse_num_items*sizeof(int));
    impulse_num = sf_read_int(impulse,impulse_buf,impulse_num_items);
    sf_close(impulse);
    printf("Read %d items\n",impulse_num);


    float *output_buf;
    int output_buf_size = num + impulse_num - 1;

    //allocate space for float data in output data;
    output_buf = (float *) malloc(output_buf_size*sizeof(float));

    //scaling items in input buffer and impulse buffer
    float * newInputBuf = scaleInputs(buf, num_items);
    float * newImpulseBuf = scaleInputs(impulse_buf, impulse_num_items);

    printf("-----------------------\n");

    printf("convolving\n");
    //convolve the current one
    convolve(newInputBuf, num, newImpulseBuf, impulse_num, output_buf, output_buf_size);
    printf("finished convolving\n");


    //int temp_output_size = 19999;
    int * final_output_buf = scaleOutputs(output_buf, output_buf_size);
    printf("here1\n");



    SNDFILE *out;
    SF_INFO outInfo;

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

    int outputSize = sf_write_int(out, final_output_buf, output_buf_size); // temp_output_size should be actual output size
    
    printf("here3\n");

    sf_close(impulse);
    sf_close(out);
    return 0;
}