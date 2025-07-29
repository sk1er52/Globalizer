/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      Calculatin_factory.h                                            //
//                                                                         //
//  Purpose:   Header file for method factory class                        //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

#ifndef __CALCULATION_FACTORY_H__
#define __CALCULATION_FACTORY_H__

#include "Calculation.h"
#include "Evolvent.h"


/// Класс необходим для правильного создания вычислителей
class CalculationFactory
{
public:
  ///Стандартный метод создания вычислителя, используется в процессе, при повторном создании возвращает уже созданный
  static Calculation* CreateCalculation(Task& _pTask, Evolvent* evolvent = 0);
  static Calculation* CreateCalculation2(Task& _pTask, Evolvent* evolvent = 0);
  
  ///Всегда создает новый вычислитель, НЕ ИСПОЛЬЗОВАТЬ В МНОГОШАГОВОЙ СХЕМЕ!
  static Calculation* CreateNewCalculation(Task& _pTask, Evolvent* evolvent = 0);
  
};

#endif
// - end of file ----------------------------------------------------------------------------------
