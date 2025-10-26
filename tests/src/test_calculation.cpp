#include "Calculation.h"
#include "OmpCalculation.h"
#include "Task.h"
#include "Problem.h"
#include "Parameters.h" // Убедитесь, что этот include есть

#include <gtest/gtest.h>
#include <vector>
#include <functional>

/**
 * @brief Этот тест проверяет только базовую возможность создать и корректно
 * удалить всю цепочку объектов, используемых в тестах вычислителей.
 * Если этот тест падает, проблема в конструкторах или деструкторах.
 */
TEST(CalculationSanityCheck, CreationAndDestruction)
{
    // --- Шаг 1: Полностью имитируем логику SetUp ---

    // Задаем глобальные параметры
    const int n = 2;
    parameters.Dimension = n;

    // Создаем данные для IProblem. Делаем их членами тестового класса,
    // чтобы они гарантированно жили дольше, чем 'problem'.
    std::vector<double> bounds(n, 0.0);
    std::vector<std::function<double(const double*)>> functions;
    functions.push_back([](const double* y) {
        if (y) return y[0] * 10.0;
        return 0.0;
        });

    // Создаем объекты в куче, как это делается в SetUp
    IProblem* problem = new ProblemFromFunctionPointers(n, bounds, bounds, functions, false, 0, nullptr);
    Task* task = new Task(problem, 1); // "Листовая" задача
    Calculation* ompCalc = new OMPCalculation(*task);

    // --- Шаг 2: Проверяем, что все указатели валидны ---
    ASSERT_NE(problem, nullptr);
    ASSERT_NE(task, nullptr);
    ASSERT_NE(ompCalc, nullptr);

    // --- Шаг 3: Полностью имитируем логику TearDown ---

    // Удаляем в обратном порядке
    delete ompCalc;
    delete task;
    delete problem;

    // Если мы дошли до сюда без падения, тест считается успешным
    SUCCEED();
}