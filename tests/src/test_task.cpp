#include "Task.h"
#include "ProblemManager.h"
#include "test_config.h"

#include <gtest/gtest.h>
#include <string>
#include <cstdlib>

/**
  Вспомогательный класс, помогающий задать начальную конфигурацию объекта класса #Task,
  которая будет использоваться в тестах
 */
class TaskTest : public ::testing::Test
{
protected:
  /// Размерность задачи
  static const int n = 5;
  /// Размерность подзадачи
  static const int freeN = 2;
  /// Число функционалов
  static const int numOfFunc = 1;
  /// Левая граница области поиска
  double A[MaxDim];
  /// Правая граница области поиска
  double B[MaxDim];
  /// Указатель на задачу
  Task* task;
  //Функционалы
  IProblem* problem;
  ProblemManager manager;

  void SetUp()
  {
    std::string libPath = std::string(TESTDATA_BIN_PATH) + std::string(LIB_RASTRIGIN);
    if (ProblemManager::OK_ == manager.LoadProblemLibrary(libPath))
    {
      problem = manager.GetProblem();
      problem->SetDimension(n);
      task = new Task(n, freeN, problem, 0);
    }
    else
      task = NULL;
  }

  void TearDown()
  {
    delete task;
  }
};

/**
 * Проверка параметра размерности задачи N
 * 1<= N <= MaxDim
 */
TEST_F(TaskTest, throws_when_create_with_negative_N)
{
  ASSERT_ANY_THROW(Task testTask(-1, freeN, problem, 0));
}

TEST_F(TaskTest, throws_when_create_with_null_N)
{
  ASSERT_ANY_THROW(Task testTask(0, freeN, problem, 0));
}

TEST_F(TaskTest, throws_when_create_with_too_large_N)
{
  ASSERT_ANY_THROW(Task testTask(MaxDim + 1, freeN, problem, 0));
}

TEST_F(TaskTest, throws_when_create_with_free_N_large_N)
{
  ASSERT_ANY_THROW(Task testTask(n, n + 1, problem, 0));
}

/**
 * Проверка параметра размерности подзадачи FreeN
 * 1<= FreeN <= N
 */

TEST_F(TaskTest, throws_when_create_with_negative_free_N)
{
  ASSERT_ANY_THROW(Task testTask(n, -1, problem, 0));
}

TEST_F(TaskTest, throws_when_create_with_null_free_N)
{
  ASSERT_ANY_THROW(Task testTask(n, 0, problem, 0));
}

/**
 * Создание задачи с корректными входными параметрами
 */
TEST_F(TaskTest, can_create_with_correct_values)
{
  ASSERT_NO_THROW(Task testTask(n, freeN, problem, 0));
}

/**
 * Проверка функции #SetFixed
 * FixedN == N - FreeN; FixedY != NULL
 */
TEST_F(TaskTest, throws_when_set_invalid_FixedN_value)
{
  ASSERT_FALSE(task == NULL);
  ASSERT_ANY_THROW(task->SetFixed(n - freeN + 1, A));
}

TEST_F(TaskTest, throws_when_set_null_FixedY)
{
  ASSERT_FALSE(task == NULL);
  ASSERT_ANY_THROW(task->SetFixed(n - freeN, NULL));
}

TEST_F(TaskTest, set_Fixed_with_correct_value)
{
  ASSERT_FALSE(task == NULL);
  ASSERT_NO_THROW(task->SetFixed(n - freeN, A));
}
