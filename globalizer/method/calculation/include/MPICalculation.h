/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      MPICalculation.h                                            //
//                                                                         //
//  Purpose:   Header file for MPI calculation class                       //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
 * \file MPICalculation.h
 *
 * \authors Лебедев И.
 * \date 2021
 * \copyright ННГУ им. Н.И. Лобачевского
 *
 * \brief Объявление класса #MPICalculation для синхронных MPI вычислений
 */

#ifndef __MPI_CALCULATION_H__
#define __MPI_CALCULATION_H__

#include "Calculation.h"

 /**
  * \brief Реализация вычислителя с использованием технологии MPI (синхронный режим).
  *
  * Позволяет выполнять распределенные вычисления на нескольких процессах (узлах).
  * В данном режиме главный процесс отправляет блок задач и ждет, пока все
  * рабочие процессы завершат их выполнение.
  */
class MPICalculation : public Calculation
{
protected:
	/// Флаг, указывающий, что это первый запуск вычислений
	bool isFirst = true;

	/**
	 * \brief Внутренний метод, запускающий синхронные MPI-вычисления.
	 * \param[in] inputSet Входные данные.
	 * \param[out] outputSet Выходные данные.
	 */
	void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

	/**
	 * \brief Выполняет вычисления на границах отрезка (особый случай).
	 * \param[in] inputSet Входные данные.
	 * \param[out] outputSet Выходные данные.
	 */
	void StartCalculateInBorder(InformationForCalculation& inputSet, TResultForCalculation& outputSet);

public:
	/**
	 * \brief Конструктор.
	 * \param _pTask Ссылка на объект задачи.
	 */
	MPICalculation(Task& _pTask) : Calculation(_pTask)
	{
	}

	/**
	 * \brief Выполняет вычисления для набора испытаний с использованием MPI.
	 * \param[in] inputSet Входные данные.
	 * \param[out] outputSet Выходные данные.
	 */
	virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
};

#endif
// - end of file ----------------------------------------------------------------------------------