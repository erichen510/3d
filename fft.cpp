#include "stdafx.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <conio.h>
#include <math.h>
#include <time.h>
#include "fft.h"						//Kayton: define the structure of cmpx

#define pi 3.1416

void twiddle();
//void FFT();
//void IFFT();
//void TWO_D_FFT();

extern int N;
extern struct cmpx P[5000],w[4096];


void twiddle()							//Complex roots of unity used as coefficients in the Fast Fourier Transform algorithm are known as twiddle factors
{
  int i;
  //double step,angle,amp;

  //step = 6.2832/N;						//step=2*pi/N   where N=H_COL*2
  //angle = 3.1416/2;						//angel=pi/2
  for(i=0;i<N/2;i++)						//N/2=H_COL		==>		There are H_COL data
  {
    //w[i].real = sin(angle);				//= sin(pi/2 + i*step)			= cos(i*step)
    //w[i].imag = sin(angle+1.5708);		//= sin(pi/2 + i*step + pi/2)	= -sin(i*step)

	w[i].real = cos(2*pi/N*i);				//Kayton: When i increase, w[i].real goes to 0
    w[i].imag = -sin(2*pi/N*i);				//Kayton

    /*printf("\nw[%d].real=%g, w[%d].img=%g",i,w[i].real,i,w[i].imag);getche();*/
    //angle += step;
  }
}


void FFT()						/*input sample array, # of points      */
  {
  struct cmpx temp1,temp2;      /*temporary storage variables          */
  int i,j,k,Nm;                 /*loop counter variables               */
  int upper_leg, lower_leg;     /*index of upper/lower butterfly leg   */
  int leg_diff;                 /*difference between upper/lower leg   */
  int num_stages=0;             /*number of FFT stages, or iterations  */
  int index, step;              /*index and step between twiddle factor*/


  //printf("\nStarting FFT");

/* log(base 2) of # of points = # of stages  */
  i=1;
  Nm = N;
  //printf("\nNm=%d",Nm);
  do
    {
    num_stages+=1;
    i=i*2;
	//printf("\ni=%d",i);getche();
    } while (i!=Nm);



/* starting difference between upper and lower butterfly legs*/
  leg_diff=Nm/2;													//Kayton: leg_diff=H_COL

/* step between values in twiddle factor array twiddle.h     */
  step=N/Nm;														//Kayton: step=1

/* For Nm-point FFT                                          */

  for (i=0;i<num_stages;i++)
  {
    index=0; 
    for (j=0;j<leg_diff;j++)
    {
      for (upper_leg=j; upper_leg<Nm; upper_leg+=(2*leg_diff))
      {
        lower_leg=upper_leg+leg_diff;

        temp1.real=(P[upper_leg]).real + (P[lower_leg]).real;
        temp1.imag=(P[upper_leg]).imag + (P[lower_leg]).imag;

        temp2.real=(P[upper_leg]).real - (P[lower_leg]).real;
        temp2.imag=(P[upper_leg]).imag - (P[lower_leg]).imag;

        (P[lower_leg]).real=temp2.real*(w[index]).real-temp2.imag*(w[index]).imag;
        (P[lower_leg]).imag=temp2.real*(w[index]).imag+temp2.imag*(w[index]).real;

        (P[upper_leg]).real=temp1.real;
        (P[upper_leg]).imag=temp1.imag;
      }
      index+=step;
    }

    leg_diff=leg_diff/2;
    step*=2;
  }

/* bit reversal for resequencing data */
  j=0;
  for (i=1;i<(Nm-1);i++)
  {
    k=Nm/2;
    while (k<=j)
    {
      j=j-k;
      k=k/2;
    }

    j=j+k;

    if (i<j)
    {
      temp1.real=(P[j]).real;
      temp1.imag=(P[j]).imag;
      (P[j]).real=(P[i]).real;
      (P[j]).imag=(P[i]).imag;
      (P[i]).real=temp1.real;
      (P[i]).imag=temp1.imag;
     }
  }
  return;
}



void IFFT() /*input sample array, # of points      */
  {
  struct cmpx temp1,temp2;      /*temporary storage variables          */
  int i,j,k,Nm;                /*loop counter variables               */
  int upper_leg, lower_leg; /*index of upper/lower butterfly leg   */
  int leg_diff;             /*difference between upper/lower leg   */
  int num_stages=0;         /*number of FFT stages, or iterations  */
  int index, step;          /*index and step between twiddle factor*/

  Nm = N;
  for(i=0;i<Nm;i++)
    P[i].imag = -P[i].imag;

/* log(base 2) of # of points = # of stages  */
  i=1;

  do
    {
    num_stages+=1;
    i=i*2;
    } while (i!=Nm);

/* starting difference between upper and lower butterfly legs*/
  leg_diff=Nm/2;
/* step between values in twiddle factor array twiddle.h     */
  step=N/Nm;
/* For Nm-point IFFT                                           */

  for (i=0;i<num_stages;i++)
    {
    index=0;
    for (j=0;j<leg_diff;j++)
      {
      for (upper_leg=j;upper_leg<Nm;upper_leg+=(2*leg_diff))
        {
        lower_leg=upper_leg+leg_diff;
        temp1.real=(P[upper_leg]).real + (P[lower_leg]).real;
		temp1.imag=(P[upper_leg]).imag + (P[lower_leg]).imag;
        temp2.real=(P[upper_leg]).real - (P[lower_leg]).real;
		temp2.imag=(P[upper_leg]).imag - (P[lower_leg]).imag;
		(P[lower_leg]).real=temp2.real*(w[index]).real-temp2.imag*(w[index]).imag;
		(P[lower_leg]).imag=temp2.real*(w[index]).imag+temp2.imag*(w[index]).real;
        (P[upper_leg]).real=temp1.real;
        (P[upper_leg]).imag=temp1.imag;
        }
      index+=step;
      }
    leg_diff=leg_diff/2;
    step*=2;
    }
/* bit reversal for resequencing data */
  j=0;
  for (i=1;i<(Nm-1);i++)
    {
    k=Nm/2;
    while (k<=j)
      {
      j=j-k;
      k=k/2;
      }
    j=j+k;
    if (i<j)
      {
      temp1.real=(P[j]).real;
      temp1.imag=(P[j]).imag;
      (P[j]).real=(P[i]).real;
      (P[j]).imag=(P[i]).imag;
      (P[i]).real=temp1.real;
      (P[i]).imag=temp1.imag;
      }
    }

  for(i=0;i<Nm;i++)
  {
    P[i].real /= Nm;
    P[i].imag /= -Nm;
  }

  return;
}

void FFT_SHIFT()
{
	int i,k,p;
	double m,n;

	k=N/2;
	for(i=0;i<k;i++)
	{
		p=i+k;
		m=P[i].real ;
		n=P[i].imag ;
		P[i].real = P[p].real;
		P[i].imag = P[p].imag;
		P[p].real = m;
		P[p].imag = n;
	}
}