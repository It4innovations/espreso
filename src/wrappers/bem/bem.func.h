// Dalibor Lukas, March 2008


#ifndef FFFUNC_H
#define FFFUNC_H

#include <math.h>
#include <complex>
#include <string.h>
#include <stdlib.h>
#include <limits.h>

#define double_complex complex<double>

extern double_complex Imag;

template<class C>
inline C mymin (C x, C y)
{
  return x<y ? x : y;
};

template<class C>
inline C mymax (C x, C y)
{
  return x>y ? x : y;
};

//template<class C>
//inline double abs (C x)
//{
//  return x>0 ? x : -x;
//};

inline double cabs (double_complex x)
{
  return sqrt(real(x)*real(x)+imag(x)-imag(x));
};

template<class C>
inline int sign (C x)
{
  return x>0 ? 1 : ( x<0 ? -1 : 0 );
};


template<class C> double myabs(C val);
template<class C> double myreal(C val);
template<class C> double myimag(C val);
template<class C> C myconj(C val);


typedef struct {double x; double y;} Point2D;
typedef struct {double x; double y; double z;} Point3D;


void strcat_itoa(char *string, unsigned int number);


template<class C>
void BubbleSort (int n, C *array, int ascending=1, int *indices=0)
{
  int i, flag;
  C tmp;
  int tmpIdx;
  if (ascending)
    do
      {
        flag = 1;
        for (i=0; i<n-1; i++)
	  if (array[i] > array[i+1])
	    {
	      tmp = array[i];
	      array[i] = array[i+1];
	      array[i+1] = tmp;
	      flag = 0;
              if (indices)
                {
                  tmpIdx = indices[i];
                  indices[i] = indices[i+1];
                  indices[i+1] = tmpIdx;
                }
	    }
      }
    while (!flag);
  else
    do
      {
        flag = 1;
        for (i=0; i<n-1; i++)
	  if (array[i] < array[i+1])
	    {
	      tmp = array[i];
	      array[i] = array[i+1];
	      array[i+1] = tmp;
	      flag = 0;
              if (indices)
                {
                  tmpIdx = indices[i];
                  indices[i] = indices[i+1];
                  indices[i+1] = tmpIdx;
                }
	    }
      }
    while (!flag);
};


template<class C>
void QuickSort (C* array, int l, int r,
                int ascending=1, int *indices=0)
{
  int i, j, tmpIdx;
  C pivot, tmp;

  i = l;
  j = r;
  pivot = array[(l+r)/2];

  do
    {
      for ( ; i<r && (ascending ? array[i]<pivot : array[i]>pivot) ; i++);
      for ( ; j>l && (ascending ? pivot<array[j] : pivot>array[j]) ; j--);
      if (i<=j)
        {
          if (i<j)
            {
              tmp = array[i];
              array[i] = array[j];
              array[j] = tmp;
              if (indices)
                {
                  tmpIdx = indices[i];
                  indices[i] = indices[j];
                  indices[j] = tmpIdx;
                }
            }
          i++;
          j--;
        }
    }
  while (i<=j);

  if (j>l)
    QuickSort<C>(array,l,j,ascending,indices);
  if (i<r)
    QuickSort<C>(array,i,r,ascending,indices);
};


// Add elem to elems and increase nelems while preserving it unique and sorted
template<class C>
int InsertSorted (C elem, C *elems, int &nelems)
{
  int i, i0, i1;
  if (nelems==0)
    {
      elems[nelems++] = elem;
      return 1;
    }
  if (elem<elems[0])
    {
      memmove(elems+1,elems,nelems*sizeof(C));
      elems[0] = elem;
      nelems++;
      return 1;
    }
  if (elem>elems[nelems-1])
    {
      elems[nelems++] = elem;
      return 1;
    }
  
  i0 = 0;
  i = nelems/2;
  i1 = nelems-1;
  while (i1-i0>1)
    {
      if (elem==elems[i])
	return 0;
      if (i==0 && elem<elems[i])
	{
	  memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
	  elems[i] = elem;
	  nelems++;
	  return 1;
	}
      if (i>0 && elems[i-1]<elem && elem<elems[i])
	{
	  memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
	  elems[i] = elem;
	  nelems++;
	  return 1;
	}      
      if (elem<elems[i])
	i1 = i;
      else
	i0 = i;
      i = i0 + (i1-i0)/2;
    }
  if (elem==elems[i])
    return 0;
  if (i==0 && elem<elems[i])
    {
      memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
      elems[i] = elem;
      nelems++;
      return 1;
    }
  if (i>0 && elems[i0]<elem && elem<elems[i0+1])
    {
      memmove(elems+i0+2,elems+i0+1,(nelems-(i0+1))*sizeof(C));
      elems[i0+1] = elem;
      nelems++;
      return 1;
    }      
  if (i1==nelems-1 && elem>elems[i1])
    {
      i = i1;
      elems[nelems++] = elem;
      return 1;
    }    
  return 0;
};	  


// Add newElems to elems and increase nelems while preserving unique and sorted
template<class C>
int InsertSorted (C *newElems, int nnewElems, C *elems, int &nelems)
{
  int i;
  for (i=0; i<nnewElems; i++)
    InsertSorted<C>(newElems[i],elems,nelems);
  return 0;
};


// Add elem to nonunique elems and increase nelems while preserving ascending sorted.
// If indices are nonzero, idx is inserted at the same position.
template<class C>
void InsertNonuniqueSorted (C elem, C *elems, int &nelems, int ascending=1,
                            int *indices=0, int idx=0)
{
  int i, i0, i1;
  if (nelems==0)
    {
      if (indices)
        indices[nelems] = idx;
      elems[nelems++] = elem;
      return;
    }
  
  i0 = 0;
  i = nelems/2;
  i1 = nelems-1;
  while (i1-i0>1)
    {
      if ( i==0 &&
           ( ascending ? elem<=elems[i] : elem>=elems[i] ) )
        {
	  memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
	  elems[i] = elem;
          if (indices)
            {
	      memmove(indices+i+1,indices+i,(nelems-i)*sizeof(int));
	      indices[i] = idx;
            }
	  nelems++;
	  return;
	}
      if ( i>0 &&
           ( ascending ? (elems[i-1]<=elem && elem<=elems[i])
                       : (elems[i-1]>=elem && elem>=elems[i]) ) )
	{
	  memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
	  elems[i] = elem;
          if (indices)
            {
	      memmove(indices+i+1,indices+i,(nelems-i)*sizeof(int));
	      indices[i] = idx;
            }
	  nelems++;
	  return;
	}      
      if ( ascending ? elem<elems[i] : elem>elems[i] )
	i1 = i;
      else
	i0 = i;
      i = i0 + (i1-i0)/2;
    }
  if ( i==0 &&
       ( ascending ? elem<=elems[i] : elem>=elems[i] ) )
    {
      memmove(elems+i+1,elems+i,(nelems-i)*sizeof(C));
      elems[i] = elem;
      if (indices)
        {
          memmove(indices+i+1,indices+i,(nelems-i)*sizeof(int));
          indices[i] = idx;
        }
      nelems++;
      return;
    }
  if ( i>0 &&
       ( ascending ? (elems[i0]<=elem && elem<=elems[i0+1])
                   : (elems[i0]>=elem && elem>=elems[i0+1]) ) )
    {
      memmove(elems+i0+2,elems+i0+1,(nelems-(i0+1))*sizeof(C));
      elems[i0+1] = elem;
      if (indices)
        {
          memmove(indices+i0+2,indices+i0+1,(nelems-(i0+1))*sizeof(int));
          indices[i0+1] = idx;
        }
      nelems++;
      return;
    }      
  if ( i1==nelems-1 &&
       ( ascending ? elem>=elems[i1] : elem<=elems[i1] ) )
    {
      i = i1;
      if (indices)
        indices[nelems] = idx;
      elems[nelems++] = elem;
      return;
    }    
};	  


// Find elem in elems and return its index (indexed from 1) or 0 otherwise
template<class C>
int Find (C elem, C *elems, int nelems)
{
  int i;
  for (i=0; i<nelems; i++)
    if (elems[i]==elem)
      return i+1;
  return 0;
};


// Insert (sorted, unique) idx2 to DOFsToDOFs1
// and increase nDOFsToDOFs1
void Insert (int *DOFsToDOFs1, int *nDOFsToDOFs1, int idx2);

void InsertSymmetric (int idx1, int *DOFsToDOFs1,
  	              int *nDOFsToDOFs1, int idx2,
		      int *DOFsToDOFs2, int *nDOFsToDOFs2);


#endif
