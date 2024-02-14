// Dalibor Lukas, March 2008


using namespace std;

#include "bem.func.h"


double_complex Imag(0,1);

template<>
double myabs<int> (int val)
{
  return (val>0) ? val : -val;
}

/*
template<>
float myabs<float> (float val)
{
  return fabs(val);
}
*/

template<>
double myabs<double> (double val)
{
  return fabs(val);
}


template<>
double myabs<double_complex> (double_complex val)
{
  return cabs(val);
}


/*
template<>
float myreal<float> (float val)
{
  return val;
}
*/

template<>
double myreal<double> (double val)
{
  return val;
}


template<>
double myreal<double_complex> (double_complex val)
{
  return real(val);
}


template<>
double myimag<double> (double val)
{
  return 0.0;
}


template<>
double myimag<double_complex> (double_complex val)
{
  return imag(val);
}

/*
template<>
float myconj<float> (float val)
{
  return val;
}
*/

template<>
double myconj<double> (double val)
{
  return val;
}


template<>
double_complex myconj<double_complex> (double_complex val)
{
  return conj(val);
}


void strcat_itoa(char *string, unsigned int number)
{
  unsigned int n = number;
  int order = 0;
  while (n/10 > 0)
    {
      n /= 10;
      order++;
    }
  string[order+1] = 0;
  while (order>=0)
    {
      string[order] = '0'+(number%10);
      number /= 10;
      order--;
    }
}


// Insert (sorted, unique) idx2 to DOFsToDOFs1
// and increase nDOFsToDOFs1
void Insert (int *DOFsToDOFs1, int *nDOFsToDOFs1, int idx2)
{
  int i, n;
  n = *nDOFsToDOFs1;
  if (n==0)
    {
      DOFsToDOFs1[0] = idx2;
      (*nDOFsToDOFs1) ++;
    }
  else if (abs(idx2)<abs(DOFsToDOFs1[0]))
    {
      memmove(DOFsToDOFs1+1,DOFsToDOFs1,n*sizeof(int));
      DOFsToDOFs1[0] = idx2;
      (*nDOFsToDOFs1) ++;
    }
  else if (abs(idx2)>abs(DOFsToDOFs1[n-1]))
    {
      DOFsToDOFs1[n] = idx2;
      (*nDOFsToDOFs1) ++;
    }
  else
    for (i=0; i<n-1; i++)
      if (abs(idx2)>abs(DOFsToDOFs1[i]) && abs(idx2)<abs(DOFsToDOFs1[i+1]))
	{
	  memmove(DOFsToDOFs1+i+2,DOFsToDOFs1+i+1,(n-(i+1))*sizeof(int));
	  DOFsToDOFs1[i+1] = idx2;
	  (*nDOFsToDOFs1) ++;
	}
}


void InsertSymmetric (int idx1, int *DOFsToDOFs1,
  	              int *nDOFsToDOFs1, int idx2,
		      int *DOFsToDOFs2, int *nDOFsToDOFs2)
{
  Insert(DOFsToDOFs1,nDOFsToDOFs1,idx2);
  Insert(DOFsToDOFs2,nDOFsToDOFs2,idx1);
}


