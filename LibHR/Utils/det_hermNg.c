/*******************************************************************************
*
* File  det_hermNg.c
*
* Function to calculate determinant of an hermitian-matrix using
* LU-decomposition. See NR.
*
* Ari Hietanen
*
*******************************************************************************/

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include "utils.h"
#include "suN.h"
#include "representation.h"

#ifndef GAUGE_SON

#ifdef WITH_QUATERNIONS
void det_hermNg(double complex *res, suNg *a)
{
  suNg b;
  double det;

  b = *a;
  _suNg_quat_det(det, b);
  *res = det;
}

#else

void det_hermNg(double complex *res, suNg *a)
{
  int indx[NG];
  double d;
  int i;
  double complex tmp;
#ifdef GAUGE_SPN
  // TODO: Simplifications may be possible for SPN
  suNgfull b;
  _suNg_expand(b,*a);
#else
  suNg b=*a;
#endif
  ludcmp(b.c, indx, &d, NG);
  *res = d;
  for (i = 0; i < NG; ++i)
  {
    _complex_mul(tmp, *res, b.c[NG * i + i]);
    *res = tmp;
  }
}

#endif //WITH_QUATERNIONS

#else

/* The incoming matrix will be destroyed */

void det_hermNg(double *res, suNg *a)
{
  double vv[NG];
  double csum, ctmp, cdum;
  double big, dum, tmp;
  int i, j, k, imax;
  int d = 1;
  for (j = 0; j < NG; j++)
  {
    big = 0.0;
    for (i = 0; i < NG; i++)
    {
      tmp = a->c[NG * j + i] * a->c[NG * j + i];
      if (tmp > big)
        big = tmp;
    }
    if (big == 0.0)
    {
      *res = 0.;
      return;
    }
    vv[j] = 1.0 / sqrt(big);
  }
  imax = 0;
  for (j = 0; j < NG; j++)
  {
    for (i = 0; i < j; i++)
    {
      csum = a->c[NG * i + j];
      for (k = 0; k < i; k++)
      {
        ctmp = a->c[NG * i + k] * a->c[NG * k + j];
        csum -= ctmp;
      }
      a->c[NG * i + j] = csum;
    }
    big = 0.0;
    for (i = j; i < NG; i++)
    {
      csum = a->c[NG * i + j];
      for (k = 0; k < j; k++)
      {
        ctmp = a->c[NG * i + k] * a->c[NG * k + j];
        csum -= ctmp;
      }
      a->c[NG * i + j] = csum;
      dum = vv[i] * fabs(csum);
      if (dum >= big)
      {
        big = dum;
        imax = i;
      }
    }

    if (j != imax)
    {
      for (k = 0; k < NG; k++)
      {
        cdum = a->c[NG * imax + k];
        a->c[NG * imax + k] = a->c[NG * j + k];
        a->c[NG * j + k] = cdum;
      }
      d = -d;
      vv[imax] = vv[j];
    }

    if (fabs(a->c[NG * j + j]) == 0.0)
    {
      a->c[NG * j + j] = 1e-20;
    }

    if (j != NG - 1)
    {
      cdum = 1. / a->c[NG * j + j];
      for (i = j + 1; i < NG; i++)
      {
        ctmp = a->c[NG * i + j] * cdum;
        a->c[NG * i + j] = ctmp;
      }
    }
  }
  csum = d;
  for (j = 0; j < NG; j++)
  {
    ctmp = a->c[NG * j + j] * csum;
    csum = ctmp;
  }
  *res = csum;
}
#endif
