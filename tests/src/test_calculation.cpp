#include "Calculation.h"
#include "OmpCalculation.h"
#include "Task.h"
#include "Problem.h"
#include "Parameters.h" // ���������, ��� ���� include ����

#include <gtest/gtest.h>
#include <vector>
#include <functional>

/**
 * @brief ���� ���� ��������� ������ ������� ����������� ������� � ���������
 * ������� ��� ������� ��������, ������������ � ������ ������������.
 * ���� ���� ���� ������, �������� � ������������� ��� ������������.
 */
TEST(CalculationSanityCheck, CreationAndDestruction)
{
    // --- ��� 1: ��������� ��������� ������ SetUp ---

    // ������ ���������� ���������
    const int n = 2;
    parameters.Dimension = n;

    // ������� ������ ��� IProblem. ������ �� ������� ��������� ������,
    // ����� ��� �������������� ���� ������, ��� 'problem'.
    std::vector<double> bounds(n, 0.0);
    std::vector<std::function<double(const double*)>> functions;
    functions.push_back([](const double* y) {
        if (y) return y[0] * 10.0;
        return 0.0;
        });

    // ������� ������� � ����, ��� ��� �������� � SetUp
    IProblem* problem = new ProblemFromFunctionPointers(n, bounds, bounds, functions, false, 0, nullptr);
    Task* task = new Task(problem, 1); // "��������" ������
    Calculation* ompCalc = new OMPCalculation(*task);

    // --- ��� 2: ���������, ��� ��� ��������� ������� ---
    ASSERT_NE(problem, nullptr);
    ASSERT_NE(task, nullptr);
    ASSERT_NE(ompCalc, nullptr);

    // --- ��� 3: ��������� ��������� ������ TearDown ---

    // ������� � �������� �������
    delete ompCalc;
    delete task;
    delete problem;

    // ���� �� ����� �� ���� ��� �������, ���� ��������� ��������
    SUCCEED();
}