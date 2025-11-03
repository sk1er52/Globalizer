#pragma once

#include "ProblemInterface.h"
#include "Parameters.h"

/**
 * \brief Класс-заглушка для интерфейса IProblem.
 *
 * \details Предоставляет предсказуемую и контролируемую реализацию интерфейса IProblem,
 * необходимую для изоляции тестируемых классов (Calculation, OMPCalculation и др.)
 * от реальных, сложных реализаций задач оптимизации.
 */

class CalculationProblem : public IProblem {
private:
    /// Внутреннее хранилище для размерности задачи.
    int m_dimension;

public:
    /**
     * \brief Конструктор.
     * \details Инициализирует мок-объект, используя глобальное значение размерности.
     */
    CalculationProblem() {
        m_dimension = parameters.Dimension;
    }

    /**
     * \brief Виртуальный деструктор.
     */
    virtual ~CalculationProblem() {}

    // --- Реализация чисто виртуальных методов для соответствия интерфейсу ---
    // Эти методы необходимы, чтобы компилятор разрешил создавать объекты CalculationProblem.
    // Их реализация минимальна, так как они не влияют на логику тестируемых классов.

    int SetConfigPath(const std::string& configPath) override {
        return IProblem::OK;
    }
     
    int SetDimension(int dimension) override {
        m_dimension = dimension;
        return IProblem::OK;
    }

    int GetDimension() const override {
        return m_dimension;
    }

    int Initialize() override {
        return IProblem::OK;
    }

    int GetNumberOfConstraints() const override {
        return 1;
    }

    int GetNumberOfCriterions() const override {
        return 1;
    }

    // --- Методы, которые активно используются в тестах ---

    /**
     * \brief Возвращает общее число функций.
     * \details В нашем CalculationProblem это всегда 2 (1 ограничение + 1 критерий).
     */
    int GetNumberOfFunctions() const override {
        return 2;
    }

    /**
     * \brief Задает границы области поиска.
     */
    void GetBounds(double* lower, double* upper) override {
        for (int i = 0; i < m_dimension; ++i) {
            lower[i] = -1.0;
            upper[i] = 1.0;
        }
    }

    /**
     * \brief Вычисляет значение функции по предсказуемой формуле.
     * \details Возвращает (сумма координат - 10) для первой функции и
     * (сумма координат + 10) для второй. Это позволяет тестам проверять
     * корректность вычислений.
     */
    double CalculateFunctionals(const double* y, int fNumber) override {
        double sum = 0.0;
        for (int i = 0; i < m_dimension; ++i) {
            sum += y[i];
        }
        return fNumber == 0 ? sum - 10.0 : sum + 10.0;
    }

    int GetOptimumValue(double& value) const override {
        return IProblem::UNDEFINED;
    }

    int GetOptimumPoint(double* y) const override {
        return IProblem::UNDEFINED;
    }
};