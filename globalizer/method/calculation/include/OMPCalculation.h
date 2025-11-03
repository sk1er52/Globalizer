/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      OMPCalculation.h                                            //
//                                                                         //
//  Purpose:   Header file for OpenMP calculation class                    //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
 * \file OMPCalculation.h
 *
 * \authors Лебедев И.
 * \date 2015
 * \copyright ННГУ им. Н.И. Лобачевского
 *
 * \brief Объявление класса #OMPCalculation для вычислений с помощью OpenMP
 */

#ifndef __OMP_CALCULATION_H__
#define __OMP_CALCULATION_H__

#include "Calculation.h"

 /**
  * \brief Реализация вычислителя с использованием технологии OpenMP.
  *
  * Позволяет распараллеливать вычисления на несколько потоков
  * в рамках одного процесса.
  */
class OMPCalculation : public Calculation
{
protected:
	/**
	 * \brief Внутренний метод, запускающий процесс вычислений.
	 * \param[in] inputSet Входные данные.
	 * \param[out] outputSet Выходные данные.
	 */
	void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

public:
	/**
	 * \brief Конструктор.
	 * \param _pTask Ссылка на объект задачи.
	 */
	OMPCalculation(Task& _pTask) : Calculation(_pTask)
	{
	}

	/**
	 * \brief Выполняет вычисления для набора испытаний с использованием OpenMP.
	 * \param[in] inputSet Входные данные.
	 * \param[out] outputSet Выходные данные.
	 */
	virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
};

#endif
// - end of file ----------------------------------------------------------------------------------