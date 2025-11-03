/////////////////////////////////////////////////////////////////////////////
//                                                                         //
//             LOBACHEVSKY STATE UNIVERSITY OF NIZHNY NOVGOROD             //
//                                                                         //
//                       Copyright (c) 2015 by UNN.                        //
//                          All Rights Reserved.                           //
//                                                                         //
//  File:      InformationForCalculation.h                                 //
//                                                                         //
//  Purpose:   Header file for method class                                //
//                                                                         //
//  Author(s):                                                             //
//                                                                         //
/////////////////////////////////////////////////////////////////////////////

/**
 * \file InformationForCalculation.h
 *
 * \authors
 * \date
 * \copyright ННГУ им. Н.И. Лобачевского
 *
 * \brief Объявление структур данных для вычислителей
 *
 * \details Содержит объявление структур #InformationForCalculation и
 * #TResultForCalculation, используемых для передачи данных в
 * вычислители и получения результатов.
 */


#ifndef __INFORMATION_FOR_CALCULATION_H__
#define __INFORMATION_FOR_CALCULATION_H__

#include <vector>
#include <fstream>
#include <algorithm>
#include "Trial.h"


 // ------------------------------------------------------------------------------------------------

 /**
  * \brief Структура для передачи входных данных в вычислитель.
  *
  * Инкапсулирует информацию, необходимую для проведения
  * одного или нескольких испытаний.
  */
class InformationForCalculation
{
public:
    /**
     * \brief Массив указателей на испытания, которые нужно вычислить.
     *
     * В зависимости от типа метода в данном массиве может быть:
     * - одна компонента (последовательный алгоритм)
     * - #NumPoints компонент (параллельный синхронный и асинхронный алгоритмы)
     */
    std::vector<Trial*> trials;

    /**
     * \brief Оператор доступа к элементам массива trials.
     * \param i Индекс элемента.
     * \return Ссылка на указатель на испытание.
     */
    Trial*& operator [] (int i)
    {
        if (i >= 0 && i < trials.size())
            return trials[i];
        else
            throw "operator [] (int i)";
    }

    /**
     * \brief Обнуляет все указатели в массиве trials.
     */
    void ToZero()
    {
        for (unsigned int i = 0; i < trials.size(); i++)
        {
            trials[i] = 0;
        }
    }

    InformationForCalculation()
    {
        trials.reserve(1);
        ToZero();
    }

    /**
     * \brief Возвращает текущий размер массива trials.
     * \return Количество испытаний.
     */
    int GetSize()
    {
        return (int)trials.size();
    }

    /**
     * \brief Изменяет размер массива trials.
     * \param n Новый размер.
     */
    void Resize(size_t n)
    {
        Clear();

        trials.resize(n);

        ToZero();
    }
    /**
     * \brief Очищает массив trials.
     */
    void Clear()
    {
        trials.clear();
    }
};

/**
 * \brief Структура для получения результатов от вычислителя.
 */
struct TResultForCalculation
{
    /**
     * \brief Массив указателей на вычисленные испытания.
     *
     * Содержит те же указатели, что и #InformationForCalculation,
     * но с уже вычисленными значениями функций.
     */
    std::vector<Trial*> trials;

    /// Количество вычислений каждой функции
    std::vector<int> countCalcTrials;

    /// Уровень процесса, на котором была вычислена точка
    std::vector<int> procLevel;

    std::vector<Trial*> NeighboursAdditionalPoints;
    std::vector<int> NeighboursAdditionalProcLevel;

    TResultForCalculation()
    {
        trials.reserve(1);
        countCalcTrials.reserve(1);
        procLevel.reserve(1);
        NeighboursAdditionalPoints.reserve(1);
        NeighboursAdditionalProcLevel.reserve(1);
    }

    /**
     * \brief Изменяет размер векторов с результатами.
     * \param n Новый размер.
     */
    void Resize(size_t n)
    {
        trials.resize(n);
        countCalcTrials.resize(n);
        procLevel.resize(n);
        NeighboursAdditionalPoints.clear();
        NeighboursAdditionalProcLevel.clear();
    }

    /**
     * \brief Очищает все векторы с результатами.
     */
    void Clear()
    {
        trials.clear();
        countCalcTrials.clear();
        procLevel.clear();
        NeighboursAdditionalPoints.clear();
        NeighboursAdditionalProcLevel.clear();
    }
};

#endif
// - end of file ----------------------------------------------------------------------------------