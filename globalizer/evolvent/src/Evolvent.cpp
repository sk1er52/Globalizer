/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      evolvent.cpp                                                //
//                                                                         //
//  Purpose:   Source file for evolvent classes                            //
//                                                                         //
//  Author(s): Barkalov K., Sysoyev A.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#pragma warning(disable:4996) 

#include <string.h>
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <iostream>

#include "Evolvent.h"
#include "Exception.h"

// ------------------------------------------------------------------------------------------------
Evolvent::Evolvent(int _N, int _m) : extNull(0.0), extOne(1.0), extHalf(0.5)
{
  int i;
  if ((_N < 1) || (_N > MaxDim))
  {
    throw EXCEPTION("N is out of range");
  }
  N = _N;
  y = new double[N];

  if ((_m < 2) || (_m > MaxM))
  {
    throw EXCEPTION("m is out of range");
  }
  m = _m;

  for (nexpExtended = extOne, i = 0; i < N; nexpExtended += nexpExtended, i++)
    ;
}

// ------------------------------------------------------------------------------------------------
Evolvent::~Evolvent()
{
  delete[] y;
}
// ------------------------------------------------------------------------------------------------
Evolvent::Evolvent(const Evolvent& evolvent) : extNull(0.0), extOne(1.0), extHalf(0.5)
{
  //Считаем развертку evolvent корректной, проверка ее параметров не требуется
  N = evolvent.N;
  y = new double[N];
  m = evolvent.m;
  nexpExtended = evolvent.nexpExtended;

  for (int i = 0; i < N; i++)
  {
    A[i] = evolvent.A[i];
    B[i] = evolvent.B[i];
  }
}

// ------------------------------------------------------------------------------------------------
Evolvent& Evolvent::operator=(const Evolvent& evolvent)
{
  //Считаем развертку evolvent корректной, проверка ее параметров не требуется
  if (N != evolvent.N)
  {
    N = evolvent.N;
    if (y)
    {
      delete[] y;
    }
    y = new double[N];
  }
  m = evolvent.m;
  nexpExtended = evolvent.nexpExtended;

  for (int i = 0; i < N; i++)
  {
    A[i] = evolvent.A[i];
    B[i] = evolvent.B[i];
  }

  return *this;
}


// ----------------------------------------------------------------------------
void Evolvent::SetBounds(const double* _A, const double* _B)
{
  for (int i = 0; i < N; i++)
  {
    A[i] = _A[i];
    B[i] = _B[i];
  }
}

// ----------------------------------------------------------------------------
//void evolvent::CalculateNumbr(Extended* s, int* u, int* v, int* l)
void Evolvent::CalculateNumbr(Extended *s, long long *u, long long *v, long long *l)
// calculate s(u)=is,l(u)=l,v(u)=iv by u=iu
{
//  int i, k1, k2, l1;
  long long i, k1, k2, l1;
  Extended is, iff;

  iff = nexpExtended;
  is = extNull;
  k1 = -1;
  k2 = 0;
  l1 = 0;
  for (i = 0; i < N; i++)
  {
    iff = iff / 2;
    k2 = -k1 * u[i];
    v[i] = u[i];
    k1 = k2;
    if (k2 < 0)
      l1 = i;
    else
    {
      is += iff;
      *l = i;
    }
  }
  if (is == extNull)
    *l = N - 1;
  else
  {
    v[N - 1] = -v[N - 1];
    if (is == (nexpExtended - extOne))
      *l = N - 1;
    else
    {
      if (l1 == (N - 1))
        v[*l] = -v[*l];
      else
        *l = l1;
    }
  }
  *s = is;
}

// ----------------------------------------------------------------------------
//void evolvent::CalculateNode(Extended is, int n, int *u, int *v, int *l)
void Evolvent::CalculateNode(Extended is, int n, long long *u, long long *v, long long *l)
// вычисление вспомогательного центра u(s) и соответствующих ему v(s) и l(s)
// calculate u=u[s], v=v[s], l=l[s] by is=s
{
//  int n1, i, j, k1, k2, iq;
  long long n1, i, j, k1, k2, iq;
  Extended iff;

  iq = 1;
  n1 = n - 1;
  *l = 0;
  if (is == 0)
  {
    *l = n1;
    for (i = 0; i < n; i++)
    {
      u[i] = -1;
      v[i] = -1;
    }
  }
  else if (is == (nexpExtended - extOne))
  {
    *l = n1;
    u[0] = 1;
    v[0] = 1;
    for (i = 1; i < n; i++)
    {
      u[i] = -1;
      v[i] = -1;
    }
    v[n1] = 1;
  }
  else
  {
    iff = nexpExtended;
    k1 = -1;
    for (i = 0; i < n; i++)
    {
      iff = iff / 2;
      if (is >= iff)
      {
        if ((is == iff) && (is != extOne))
        {
          *l = i;
          iq = -1;
        }
        is -= iff;
        k2 = 1;
      }
      else
      {
        k2 = -1;
        if ((is == (iff - extOne)) && (is != extNull))
        {
          *l = i;
          iq = 1;
        }
      }
      j = -k1 * k2;
      v[i] = j;
      u[i] = j;
      k1 = k2;
    }
    v[*l] = v[*l] * iq;
    v[n1] = -v[n1];
  }
}

// ------------------------------------------------------------------------------------------------
void Evolvent::transform_P_to_D()
{
  //if (N == 1) return;
  // transformation from hypercube P to hyperinterval D
  for (int i = 0; i < N; i++)
    y[i] = y[i] * (B[i] - A[i]) + (A[i] + B[i]) / 2;
}

// ----------------------------------------------------------------------------
void Evolvent::transform_D_to_P()
{
  //if (N == 1) return;
  // transformation from hyperinterval D to hypercube P
  for (int i = 0; i < N; i++)
    y[i] = (y[i] - (A[i] + B[i]) / 2) / (B[i] - A[i]);
}

// ----------------------------------------------------------------------------
double* Evolvent::GetYOnX(const Extended& _x)
{
  if (N == 1)
  {
    y[0] = _x.toDouble() - 0.5;
    return y;
  }

  //int iu[MaxDim];
  //int iv[MaxDim];
  long long iu[MaxDim]{};
  long long iv[MaxDim]{};
  //int l;
  long long l;
  Extended d;
  int mn;
  double r;
//  int iw[MaxDim];
  long long iw[MaxDim]{};
  //int it, i, j;
  long long it, i, j;
  Extended is;

  d = _x;
  r = 0.5;
  it = 0;
  mn = m * N;
  for (i = 0; i < N; i++)
  {
    iw[i] = 1;
    y[i] = 0.0;
  }
  for (j = 0; j < m; j++)
  {
    if (_x == extOne)
    {
      is = nexpExtended - extOne;
      d = extNull; // d = 0.0;
    }
    else
    {
      //Код из старой версии - уточнить работоспособность при N > 32
      d *= nexpExtended;
      //is = (int)d.toDouble();
      is = (long long)d.toDouble();
      d -= is;
    }
    CalculateNode(is, N, iu, iv, &l);
    i = iu[0];
    iu[0] = iu[it];
    iu[it] = i;
    i = iv[0];
    iv[0] = iv[it];
    iv[it] = i;
    if (l == 0)
      l = it;
    else if (l == it)
      l = 0;
    r *= 0.5;
    it = l;
    for (i = 0; i < N; i++)
    {
      iu[i] *= iw[i];
      iw[i] *= -iv[i];
      y[i] += r * iu[i];
    }
  }
  return y;
}

//-----------------------------------------------------------------------------
Extended Evolvent::GetXOnY()
{
  //int u[MaxDim], v[MaxDim];
  long long u[MaxDim]{}, v[MaxDim]{};
  Extended x, r1;
  if (N == 1)
  {
    x = y[0] + 0.5;
    return x;
  }

  double  r;
  //int w[MaxDim];
  long long w[MaxDim]{};
  //int l;
  long long l;
  //int i, j, it;
  long long i, j, it;
  Extended is;

  for (i = 0; i < N; i++)
    w[i] = 1;
  r = 0.5;
  r1 = extOne;
  x = extNull;
  it = 0;
  for (j = 0; j < m; j++)
  {
    r *= 0.5;
    for (i = 0; i < N; i++)
    {
      u[i] = (y[i] < 0) ? -1 : 1;
      y[i] -= r * u[i];
      u[i] *= w[i];
    }
    i = u[0];
    u[0] = u[it];
    u[it] = i;
    CalculateNumbr(&is, u, v, &l);
    i = v[0];
    v[0] = v[it];
    v[it] = i;
    for (i = 0; i < N; i++)
      w[i] *= -v[i];
    if (l == 0)
      l = it;
    else
      if (l == it)
        l = 0;
    it = l;
    r1 = r1 / nexpExtended;
    x += r1 * is;
  }
  return x;
}

//----------------------------------------------------------------------------
void Evolvent::GetImage(const Extended& x, double* _y, int EvolventNum)
{

  // в одиночной развертке evolventNum не используется
  //   введен, чтобы работал полиморфизм в множественных развертках


  if ((x.toDouble() < 0) || (x.toDouble() > 1))
  {
    throw EXCEPTION("x is out of range");
  }
  // x ---> y
  GetYOnX(x); // it saves return value to y, so no need to call operator= again

  transform_P_to_D();

  memcpy(_y, y, N * sizeof(double));
}

void Evolvent::GetInverseImage(double* _y, Extended& x)
{
  // y ---> x
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  x = GetXOnY();
}

//----------------------------------------------------------------------------
void Evolvent::GetPreimages(double* _y, Extended* x)
{
  // y ---> x
  memcpy(y, _y, N * sizeof(double));
  transform_D_to_P();
  x[0] = GetXOnY();
}

// ------------------------------------------------------------------------------------------------
/// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
double Evolvent::ZeroConstraintCalc(const double* _y, int EvolventNum)
{
  return -1;
}

// - end of file ----------------------------------------------------------------------------------
