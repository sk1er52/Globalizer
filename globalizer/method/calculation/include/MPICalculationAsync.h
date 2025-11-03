/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2021 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      MPICalculationAsync.h                                       //
//                                                                         //
//  Purpose:   Header file for Async MPI calculation class                 //
//                                                                         //
//  Author(s): Lebedev I.                                                  //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
 * \file MPICalculationAsync.h
 *
 * \authors Лебедев И.
 * \date 2021
 * \copyright ННГУ им. Н.И. Лобачевского
 *
 * \brief Объявление класса #MPICalculationAsync для асинхронных MPI вычислений
 */

#ifndef __MPI_CALCULATION_ASYNC_H__
#define __MPI_CALCULATION_ASYNC_H__

#include "Calculation.h"

 /**
  * \brief Реализация вычислителя с использованием технологии MPI (асинхронный режим).
  *
  * Главный процесс отправляет по одной задаче каждому рабочему процессу.
  * Как только какой-либо рабочий завершает задачу, он отправляет результат,
  * и главный процесс немедленно отправляет ему новую задачу, не дожидаясь остальных.
  */
class MPICalculationAsync : public Calculation
{
protected:
    /// Флаг, указывающий, что это первый запуск вычислений.
    bool isFirst = true;

    /// Вектор для хранения указателей на отправленные на вычисление испытания.
    std::vector<Trial*> vecTrials;

    /// MPI-номер потомка, закончившего решение выделенной задачи
    int ChildNumRecv;

    /// "Внутренний" номер потомка среди всех потомков данного процесса
    int ChildNum;

    /**
     * \brief Первоначальная рассылка задач всем рабочим процессам.
     */
    void FirstStartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
    /**
     * \brief Отправка новой задачи освободившемуся процессу.
     */
    void StartCalculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
    /**
     * \brief Ожидание и получение результата от любого из рабочих процессов.
     */
    void RecieveCalculatedFunctional();

public:
    /**
     * \brief Статический метод для корректного завершения работы.
     *
     * Собирает оставшиеся в работе результаты перед выходом.
     */
    static void AsyncFinilize();
    /**
     * \brief Конструктор.
     * \param _pTask Ссылка на объект задачи.
     */
    MPICalculationAsync(Task& _pTask) : Calculation(_pTask)
    {
    }

    /**
     * \brief Выполняет вычисления для набора испытаний в асинхронном режиме.
     * \param[in] inputSet Входные данные.
     * \param[out] outputSet Выходные данные.
     */
    virtual void Calculate(InformationForCalculation& inputSet, TResultForCalculation& outputSet);
};

#endif
// - end of file ----------------------------------------------------------------------------------