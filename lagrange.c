#include "common.h"
#include "lagrange.h"
#include "BasOP.h"

/******************************************************************************
*                            模块内部常量                                *
******************************************************************************/

/******************************************************************************
*                            模块内部宏定义                              *
******************************************************************************/
// when N<=20, alias become samller.
// when N>=20, noise become bigger.
#define N 3
/******************************************************************************
*                            模块内部数据类型                            *
******************************************************************************/

/******************************************************************************
*                            模块内部数据结构                            *
******************************************************************************/
struct tagLAGObject{
	int SrcLen;
	int DstLen;
	float *coe;
};
/******************************************************************************
*                            模块全局变量                                *
******************************************************************************/

/******************************************************************************
*                            模块内部变量                                *
******************************************************************************/

/******************************************************************************
*                            局部函数原型                                *
******************************************************************************/

/******************************************************************************
*                            全局函数实现                                *
******************************************************************************/
// lagrange based SRC
LAGHandle LagrangeCreate(int SrcLen, int DstLen)
{
	LAGHandle filter;
	int i, j, k;
	float *p;
	float x;
	int startX;
	filter=(LAGHandle)calloc(1,sizeof(LAGHandle));
    if (!filter)
	{
		return NULL;
	}
	filter->coe=(float *)calloc(1,DstLen * (N+1) * sizeof(float)); // coe+startX
    if (!filter->coe)
	{
		LagrangeDelete(filter);
		return NULL;
	}
	filter->SrcLen=SrcLen;
	filter->DstLen=DstLen;

	p=filter->coe;
	for(i=0; i<DstLen; i++)
	{
		x=1.0*i*(SrcLen-1)/(DstLen-1);
		startX=(int)x;
		while((startX+2)>=SrcLen)
		{
			startX-=1;
		}
		p[N]=startX;
		x=x-startX;
		for(k=0; k<N; k++)
		{
			p[k]=1.0;
			for(j=0; j<N; j++)
			{
				if(j==k)
				{
					continue;
				}
				p[k] *= (x-j)/(k-j);
			}
			if (!IsNormal(p[k]))
			{
				PRINTF("need check:%f\n", p[k]); // order(e.g.479) of polynomial is too high
			}
		}
		p += (N+1);
	}

	return filter;
}


void LagrangeDelete(LAGHandle filter)
{
    if (!filter)
	{
		return ;
	}
    if (!filter->coe)
	{
		free(filter->coe);
	}
	free(filter);
}

int LagrangeProcess(LAGHandle filter, int *inBuf, int *outBuf)
{
	int i,j;
	float *p;
	float temp;
	int startX;

    if (!filter || !inBuf || !outBuf)
	{
		return -1;
	}

	p=filter->coe;
	for(i=0; i<filter->DstLen; i++)
	{
		startX=(int)(p[N]);
		temp=0.0;
		for(j=0; j<N; j++)
		{
			temp += p[j] * inBuf[2*(startX+j)+0];
		}
		outBuf[2*i+0] = CLIP(temp, MIN_32, MAX_32);
		temp=0.0;
		for(j=0; j<N; j++)
		{
			temp += p[j] * inBuf[2*(startX+j)+1];
		}
		outBuf[2*i+1] = CLIP(temp, MIN_32, MAX_32);
		p += (N+1);
	}
	return 0;
}

LAGHandle LagrangeCreate1(int SrcLen, int DstLen)
{
	LAGHandle filter;
	int i, j, k;
	float *p;
	float x;
	filter=(LAGHandle)calloc(1,sizeof(LAGHandle));
    if (!filter)
	{
		return NULL;
	}
	filter->coe=(float *)calloc(1,DstLen * SrcLen * sizeof(float));
    if (!filter->coe)
	{
		LagrangeDelete(filter);
		return NULL;
	}
	filter->SrcLen=SrcLen;
	filter->DstLen=DstLen;

	p=filter->coe;
	for(i=0; i<DstLen; i++)
	{
		for(k=0; k<SrcLen; k++)
		{
			p[k]=1.0;
			x=1.0*i*(SrcLen-1)/(DstLen-1);
			for(j=0; j<SrcLen; j++)
			{
				if(j==k)
				{
					continue;
				}
				p[k] *= (x-j)/(k-j);
				if (i==1 &&k==27)
				{
					PRINTF("need check:%f\n", p[k]);
				}
			}
			// order(e.g.479) of polynomial is too high,lead to coe overflow and process() cost too much time.
			if (!IsNormal(p[k]))
			{
				PRINTF("need check:%f\n", p[k]);
			}
		}
		p += SrcLen;
	}

	return filter;
}


void LagrangeDelete1(LAGHandle filter)
{
    if (!filter)
	{
		return ;
	}
    if (!filter->coe)
	{
		free(filter->coe);
	}
	free(filter);
}


int LagrangeProcess1(LAGHandle filter, short int *inBuf, short int *outBuf)
{
	int i,j;
	float *p;
	float temp;

    if (!filter || !inBuf || !outBuf)
	{
		return -1;
	}

	p=filter->coe;
	for(i=0; i<filter->DstLen; i++)
	{
		temp=0.0;
		for(j=0; j<filter->SrcLen; j++)
		{
			temp += p[j] * inBuf[2*j+0];
		}
		outBuf[2*i+0] = CLIP(temp, MIN_16, MAX_16);
		temp=0.0;
		for(j=0; j<filter->SrcLen; j++)
		{
			temp += p[j] * inBuf[2*j+1];
		}
		outBuf[2*i+1] = CLIP(temp, MIN_16, MAX_16);
		p += filter->SrcLen;
	}
	return 0;
}

#include "wave.h"
void lagrange_test()
{
    WavFileHandle hWavIn;
    WavFileHandle hWavOut;
    T_FORMATCHUNK tFormatchunkIn;
    T_DATACHUNK   tDatachunkIn;
    T_FORMATCHUNK tFormatchunkOut;
    T_DATACHUNK   tDatachunkOut;
    int numSamplesPerChIn, numSamplesPerChOut, uiSamplesPerFrameIn, uiSamplesPerFrameOut, numFrames;
    char *pDataIn, *pDataOut;
	int i;
    U32 uiStartTime, uiEndTime;
	tProcStat st={0,0,~0,0, 0,0, 0,0};
	LAGHandle hLAG;

    hWavIn=WavCreate(NULL, NULL);
    if(hWavIn == NULL)
    {
        printf("WavCreate error!\n");
        return ;
    }
    WavRead(hWavIn, "../sweep2ch48k32b-3dB.wav");
    WavGetFmt(hWavIn, &tFormatchunkIn);
    WavGetData(hWavIn, &tDatachunkIn);
    tFormatchunkOut=tFormatchunkIn;
    WavSetFmt(&tFormatchunkOut, 1, 2, 44100, 32);
    WavGetSamplesPerCh(hWavIn, &numSamplesPerChIn);
#if 1 // frame process
    uiSamplesPerFrameIn = tFormatchunkIn.uiSamplesPerSec / 100; // 10ms
    uiSamplesPerFrameOut = tFormatchunkOut.uiSamplesPerSec / 100; // 10ms
#else // sample process
    uiSamplesPerFrameIn = numSamplesPerChIn;
    uiSamplesPerFrameOut = numSamplesPerChIn;
#endif
    numFrames = numSamplesPerChIn / uiSamplesPerFrameIn;
    numSamplesPerChIn = numFrames * uiSamplesPerFrameIn;
    numSamplesPerChOut = numFrames * uiSamplesPerFrameOut;
    WavSetData(&tDatachunkOut, &tFormatchunkOut, numSamplesPerChOut);
    hWavOut=WavCreate(&tFormatchunkOut, &tDatachunkOut);
    if(hWavOut == NULL)
    {
        printf("WavCreate error!\n");
        return ;
    }
    pDataIn = WavGetBuf(hWavIn);
    pDataOut = WavGetBuf(hWavOut);

	// Create
	hLAG = LagrangeCreate(uiSamplesPerFrameIn,uiSamplesPerFrameOut);
	// Process
	uiStartTime = time_ms();
	for(i=0; i<numFrames; i++)
	{
		StatEnter(&st);
		LagrangeProcess(hLAG, (I32 *)pDataIn, (I32 *)pDataOut);
		StatLeave(&st);
		pDataIn += uiSamplesPerFrameIn * tFormatchunkIn.uiBlockAlign;
		pDataOut += uiSamplesPerFrameOut * tFormatchunkOut.uiBlockAlign;
	}
	uiEndTime = time_ms();
	logprintf((LEVEL(0), "%d numFrames, Process Time=%dms\n", numFrames, uiEndTime-uiStartTime));

    // write file
    WavWrite(hWavOut, "../sweep2ch48k32b_out.wav");
    WavDelete(hWavIn);
    WavDelete(hWavOut);

    return ;
}

void lagrange_test1()
{
    WavFileHandle hWavIn;
    WavFileHandle hWavOut;
    T_FORMATCHUNK tFormatchunkIn;
    T_DATACHUNK   tDatachunkIn;
    T_FORMATCHUNK tFormatchunkOut;
    T_DATACHUNK   tDatachunkOut;
    int numSamplesPerChIn, numSamplesPerChOut, uiSamplesPerFrame, numFrames;
    char *pDataIn, *pDataOut;
	int i;
    U32 uiStartTime, uiEndTime;
	tProcStat st={0,0,~0,0, 0,0, 0,0};
	LAGHandle hLAG;

    hWavIn=WavCreate(NULL, NULL);
    if(hWavIn == NULL)
    {
        printf("WavCreate error!\n");
        return ;
    }
    WavRead(hWavIn, "../LR2ch48k16b.wav");
    WavGetFmt(hWavIn, &tFormatchunkIn);
    WavGetData(hWavIn, &tDatachunkIn);
    WavGetSamplesPerCh(hWavIn, &numSamplesPerChIn);
#if 1 // frame process
    uiSamplesPerFrame = tFormatchunkIn.uiSamplesPerSec / 100; // 10ms
#else // sample process
    uiSamplesPerFrame = numSamplesPerChIn;
#endif
    numFrames = numSamplesPerChIn / uiSamplesPerFrame;
    numSamplesPerChIn = numFrames * uiSamplesPerFrame;
    numSamplesPerChOut = numFrames * 441;

    tFormatchunkOut=tFormatchunkIn;
    WavSetFmt(&tFormatchunkOut, 1, 2, 44100, 16);
    WavSetData(&tDatachunkOut, &tFormatchunkOut, numSamplesPerChOut);
    hWavOut=WavCreate(&tFormatchunkOut, &tDatachunkOut);
    if(hWavOut == NULL)
    {
        printf("WavCreate error!\n");
        return ;
    }
    pDataIn = WavGetBuf(hWavIn);
    pDataOut = WavGetBuf(hWavOut);

	// Create
	hLAG = LagrangeCreate1(480,441);
	// Process
	uiStartTime = time_ms();
	for(i=0; i<numFrames; i++)
	{
		StatEnter(&st);
		LagrangeProcess1(hLAG, (I16 *)pDataIn, (I16 *)pDataOut);
		StatLeave(&st);
		pDataIn += uiSamplesPerFrame * tFormatchunkIn.uiBlockAlign;
		pDataOut += uiSamplesPerFrame * tFormatchunkOut.uiBlockAlign;
	}
	uiEndTime = time_ms();
	logprintf((LEVEL(0), "%d numFrames, Process Time=%dms\n", numFrames, uiEndTime-uiStartTime));

    // write file
    WavWrite(hWavOut, "../LR2ch44k16b_out.wav");
    WavDelete(hWavIn);
    WavDelete(hWavOut);

    return ;
}

/******************************************************************************
*                            局部函数实现                                *
******************************************************************************/

