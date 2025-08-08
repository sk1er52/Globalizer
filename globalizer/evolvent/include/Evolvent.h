/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      evolvent.h                                                  //
//                                                                         //
//  Purpose:   Header file for evolvent classes                            //
//                                                                         //
//  Author(s): Barkalov K., Sysoyev A.                                     //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __EVOLVENT_H__
#define __EVOLVENT_H__

#include "Common.h"
#include "Extended.h"
#include <vector>

// ------------------------------------------------------------------------------------------------
class Evolvent
{
protected:
  int      m;             // accuracy of decomposition of hypercube
  int      N;             // dimension
  double   A[MaxDim];     // left and
  double   B[MaxDim];     // right bounds of search area

  double*  y;             // y point from hypercube [-1/2, 1/2]^N

  const Extended extNull; // = 0.0; //Extended(0.0);
  const Extended extOne;  // = 1.0; //Extended(1.0);
  const Extended extHalf; // = 0.5; //Extended(0.5);
  Extended nexpExtended;

//  void CalculateNumbr(Extended* s, int* u, int* v, int* l);
  void CalculateNumbr(Extended *s, long long *u, long long *v, long long *l);
//  void CalculateNode(Extended is, int n, int* u, int* v, int* l);
  void CalculateNode(Extended is, int n, long long *u, long long *v, long long *l);
  void transform_P_to_D(); // transformation from hypercube P to hyperinterval D
  void transform_D_to_P(); // transformation from hyperinterval D to hypercube P
  double* GetYOnX(const Extended& _x);
  Extended GetXOnY();
public:
  Evolvent(int _N = 2, int _m = 10);
  Evolvent(const Evolvent& evolvent);
  virtual ~Evolvent();
  void SetBounds(const double* _A, const double* _B);
  //x-->y
  virtual void GetImage(const Extended& x, double* _y, int EvolventNum = 0);
  //y-->x
  void GetInverseImage(double* _y, Extended& x);
  //y-->x
  virtual void GetPreimages(double* _y, Extended* x);
  Evolvent& operator=(const Evolvent& evolvent);

  /// Вычисляет функцию существования точки в развертки EvolventNum для y, <0 - существует
  virtual double ZeroConstraintCalc(const double* _y, int EvolventNum = 0);

};

#endif
// - end of file ----------------------------------------------------------------------------------
