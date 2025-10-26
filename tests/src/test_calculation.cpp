// test_calculation.cpp  ── скелет
#include "CalculationFactory.h"
#include "OmpCalculation.h"
#include "Task.h"
#include "Problem.h"
#include "test_config.h"

#include <gtest/gtest.h>

class CalculationFixture : public ::testing::Test {
protected:
    static const int n = 5;
    Task* task;
    IProblem* problem;

    void SetUp() override {
        parameters.Dimension = n;
        problem = new ProblemFromFunctionPointers(
            n,
            std::vector<double>(n, -2.0),
            std::vector<double>(n, 2.0),
            { [](const double* x) {
                double s = 0; for (int i = 0; i < parameters.Dimension; i++) s += x[i] * x[i];
                return s;
            } },
            true, 0.0,
            std::vector<double>(n, 0).data()
        );
        task = new Task(problem, 1);   // ProcLevel=1 → «лист»
    }

    void TearDown() override {
        delete task;
        delete problem;
        // обнуляем синглтоны между тестами
        Calculation::leafCalculation = nullptr;
        Calculation::firstCalculation = nullptr;
    }
};

// ---------- CalculationFactory ----------
TEST_F(CalculationFixture, CreatesOMPCalculationForLeaf) {
    Calculation* calc = CalculationFactory::CreateCalculation2(*task);

    ASSERT_NE(calc, nullptr);
    EXPECT_NE(dynamic_cast<OMPCalculation*>(calc), nullptr);
}

TEST_F(CalculationFixture, ReturnsSameInstanceOnSecondCall) {
    // TODO: аналогично первому тесту, но вызовите фабрику дважды
    //       и сравните указатели (EXPECT_EQ)
}

// ---------- Calculation (статические поля) ----------
TEST_F(CalculationFixture, SetCountCalculationChangesGlobal) {
    // TODO: 1) создайте OMPCalculation calc(*task);
    //       2) вызовите calc.SetCountCalculation(5);
    //       3) проверьте, что static Calculation::countCalculation == 5;
}

// ---------- OMPCalculation::Calculate ----------
TEST_F(CalculationFixture, CalculateFillsIndexAndCounters) {
    // Подготовим один Trial
    Trial tr;
    for (int i = 0; i < n; i++) tr.y[i] = 0;           // точка в (0,0,…)
    InformationForCalculation in;
    in.trials.push_back(&tr);
    TResultForCalculation out;

    OMPCalculation calc(*task);
    // TODO: вызовите calc.Calculate(in, out);
    //       затем проверьте:
    //       1) tr.index == 0
    //       2) out.countCalcTrials.size() == 1
    //       3) out.countCalcTrials[0] == 1
}