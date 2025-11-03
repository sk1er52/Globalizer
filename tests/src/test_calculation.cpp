#include <gtest/gtest.h>

#include "Calculation.h"
#include "OMPCalculation.h"
#include "CalculationFactory.h"
#include "Task.h"
#include "Parameters.h"
#include "Trial.h"
#include "InformationForCalculation.h"
#include "CalculationProblem.h"

/**
 * \brief Тестовый класс для классов Calculation, OMPCalculation и CalculationFactory.
 *
 * \details Отвечает за создание и настройку общего тестового окружения перед каждым тестом
 * и его очистку после. Это гарантирует изоляцию тестов друг от друга.
 */
class CalculationTest : public ::testing::Test {
protected:
    /// Указатель на CalculationProblem объект задачи.
    CalculationProblem* calculationProblem;
    /// Указатель на задачу для корневого процесса (не-листовая).
    Task* nonLeafTask;
    /// Указатель на задачу для дочернего процесса (листовая).
    Task* leafTask;
    /// Хранилище для исходного состояния глобальных параметров.
    Parameters originalParameters;

    /**
     * \brief Метод настройки, вызываемый перед каждым тестом.
     * \details Инициализирует глобальные параметры, создает CalculationProblem объекты и задачи,
     * а также сбрасывает статические поля класса Calculation для обеспечения
     * чистого состояния для каждого теста.
     */
    void SetUp() override {
        int argc = 1;
        char* argv[] = { const_cast<char*>("tests"), nullptr };
        parameters.Init(argc, argv, false);

        originalParameters = parameters;
        parameters.Dimension = 2;

        calculationProblem = new CalculationProblem();
        nonLeafTask = new Task(calculationProblem, 0);
        leafTask = new Task(calculationProblem, 1);

        if (Calculation::leafCalculation) {
            delete Calculation::leafCalculation;
            Calculation::leafCalculation = nullptr;
        }
        Calculation::firstCalculation = nullptr;

        SetCountCalculation(0);
        SetStartComputingAway(true);
        GetInputCalculation().Clear();
        GetResultCalculation().Clear();
    }

    /**
     * \brief Метод очистки, вызываемый после каждого теста.
     * \details Освобождает всю выделенную в SetUp память и сбрасывает статические
     * указатели, чтобы предотвратить утечки памяти и влияние на последующие тесты.
     */
    void TearDown() override {
        delete nonLeafTask;
        delete leafTask;
        delete calculationProblem;

        GetInputCalculation().Clear();
        GetResultCalculation().Clear();

        if (Calculation::leafCalculation) {
            delete Calculation::leafCalculation;
            Calculation::leafCalculation = nullptr;
        }
        Calculation::firstCalculation = nullptr;
        parameters = originalParameters;
    }

    // --- Геттеры для доступа к protected-членам Calculation ---
    void SetStartComputingAway(bool value) {
        Calculation::isStartComputingAway = value;
    }

    void SetCountCalculation(int count) {
        Calculation::countCalculation = count;
    }

    InformationForCalculation& GetInputCalculation() {
        return Calculation::inputCalculation;
    }

    TResultForCalculation& GetResultCalculation() {
        return Calculation::resultCalculation;
    }

    int GetCountCalculation() {
        return Calculation::countCalculation;
    }

    bool IsStartComputingAway() {
        return Calculation::isStartComputingAway;
    }
};

// --- Тесты для CalculationFactory ---

/**
 * \brief Проверяет, что CreateNewCalculation корректно создает OMPCalculation для листовой задачи.
 */
TEST_F(CalculationTest, Factory_CreateNewCalculation_CreatesOMP_ForLeafTask) {
    parameters.TypeCalculation = OMP;
    Calculation* calc = CalculationFactory::CreateNewCalculation(*leafTask);

    ASSERT_NE(calc, nullptr);
    EXPECT_NE(dynamic_cast<OMPCalculation*>(calc), nullptr);

    delete calc;
}

/**
 * \brief Проверяет, что CreateNewCalculation возвращает nullptr для не-листовой задачи.
 */
TEST_F(CalculationTest, Factory_CreateNewCalculation_ReturnsNull_ForNonLeafTask) {
    parameters.TypeCalculation = OMP;
    Calculation* calc = CalculationFactory::CreateNewCalculation(*nonLeafTask);
    EXPECT_EQ(calc, nullptr);
}

/**
 * \brief Проверяет логику синглтона в методе CreateCalculation2.
 * \details Убеждается, что при повторном вызове для листовой задачи возвращается
 * тот же самый экземпляр вычислителя.
 */
TEST_F(CalculationTest, Factory_CreateCalculation2_SingletonLogic) {
    parameters.TypeCalculation = OMP;

    Calculation* calc1 = CalculationFactory::CreateCalculation2(*leafTask);
    ASSERT_NE(calc1, nullptr);
    EXPECT_NE(dynamic_cast<OMPCalculation*>(calc1), nullptr);

    Calculation* calc2 = CalculationFactory::CreateCalculation2(*leafTask);
    EXPECT_EQ(calc1, calc2);
}

/**
 * \brief Проверяет, что CreateCalculation возвращает nullptr для листовой задачи.
 */
TEST_F(CalculationTest, Factory_CreateCalculation_ReturnsNull_ForLeafTask) {
    parameters.TypeCalculation = OMP;
    Calculation* calc = CalculationFactory::CreateCalculation(*leafTask);
    EXPECT_EQ(calc, nullptr);
}

// --- Тесты для OMPCalculation ---

/**
 * \brief Тестирует стандартный режим немедленного выполнения вычислений.
 * \details Проверяет, что Calculate сразу вычисляет значения для переданного испытания
 * и корректно заполняет поля index, FuncValues и счетчики.
 */
TEST_F(CalculationTest, OMPCalculation_Calculate_ImmediateExecution) {
    OMPCalculation calc(*leafTask);
    InformationForCalculation inputSet;
    TResultForCalculation outputSet;

    Trial trial;
    trial.y[0] = 1.0;
    trial.y[1] = 2.0;
    inputSet.trials.push_back(&trial);

    SetStartComputingAway(true);
    calc.Calculate(inputSet, outputSet);

    ASSERT_EQ(outputSet.trials.size(), 1);
    Trial* resultTrial = outputSet.trials[0];

    EXPECT_EQ(resultTrial->index, 1);
    EXPECT_DOUBLE_EQ(resultTrial->FuncValues[0], -7.0);
    EXPECT_DOUBLE_EQ(resultTrial->FuncValues[1], 13.0);
    ASSERT_EQ(outputSet.countCalcTrials.size(), 2);
    EXPECT_EQ(outputSet.countCalcTrials[0], 1);
    EXPECT_EQ(outputSet.countCalcTrials[1], 1);
}

/**
 * \brief Тестирует режим отложенных вычислений.
 * \details Проверяет, что Calculate накапливает испытания при нескольких вызовах
 * и выполняет их все вместе, когда счетчик достигает нуля.
 */
TEST_F(CalculationTest, OMPCalculation_Calculate_AccumulatedExecution) {
    OMPCalculation calc(*leafTask);
    Calculation::firstCalculation = &calc;

    InformationForCalculation inputSet1, inputSet2;
    TResultForCalculation outputSet1, outputSet2;

    Trial* trial1 = new Trial();
    trial1->y[0] = 1.0;
    trial1->y[1] = 1.0;
    inputSet1.trials.push_back(trial1);

    Trial* trial2 = new Trial();
    trial2->y[0] = 2.0;
    trial2->y[1] = 2.0;
    inputSet2.trials.push_back(trial2);

    calc.SetCountCalculation(2);

    calc.Calculate(inputSet1, outputSet1);

    ASSERT_EQ(GetInputCalculation().trials.size(), 1);
    EXPECT_EQ(trial1->index, -2);

    calc.Calculate(inputSet2, outputSet2);

    ASSERT_EQ(GetResultCalculation().trials.size(), 2);
    Trial* resultTrial1 = GetResultCalculation().trials[0];
    Trial* resultTrial2 = GetResultCalculation().trials[1];

    EXPECT_EQ(resultTrial1->index, 1);
    EXPECT_DOUBLE_EQ(resultTrial1->FuncValues[0], 2.0 - 10.0);
    EXPECT_DOUBLE_EQ(resultTrial1->FuncValues[1], 2.0 + 10.0);

    EXPECT_EQ(resultTrial2->index, 1);
    EXPECT_DOUBLE_EQ(resultTrial2->FuncValues[0], 4.0 - 10.0);
    EXPECT_DOUBLE_EQ(resultTrial2->FuncValues[1], 4.0 + 10.0);

    delete trial1;
    delete trial2;
}

// --- Тесты для базового класса Calculation ---

/**
 * \brief Проверяет корректность работы сеттеров базового класса Calculation.
 */
TEST_F(CalculationTest, BaseCalculation_Setters) {
    OMPCalculation calc(*leafTask);

    Task* newTask = new Task(calculationProblem, 2);
    calc.SetTask(newTask);

    SearchData searchData(2);
    calc.SetSearchData(&searchData);

    calc.SetCountCalculation(5);
    EXPECT_EQ(GetCountCalculation(), 5);
    EXPECT_FALSE(IsStartComputingAway());

    delete newTask;
}