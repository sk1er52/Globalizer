/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      method_factory.cpp                                          //
//                                                                         //
//  Purpose:   Source file for method factory class                        //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////
#include "MethodFactory.h"

#include "CalculationFactory.h"


// ------------------------------------------------------------------------------------------------
IMethod* MethodFactory::CreateMethod(Task& _pTask, SearchData& _pData,
  Calculation& _Calculation, Evolvent& _Evolvent)
{
  IMethod* pMethod = 0;
  if ((parameters.TypeMethod == StandartMethod) ||
    (parameters.TypeMethod == HybridMethod) ||
    (parameters.TypeMethod == ManyNumPointMethod))
    pMethod = new Method(_pTask, _pData, _Calculation, _Evolvent);

  return pMethod;
}