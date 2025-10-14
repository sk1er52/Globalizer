<p align="center">
  <img src="/docs/iOpt_logo.png" width="200" height="150"/>
</p>

[![License: BSD 3-Clause](https://img.shields.io/badge/License-BSD%203--Clause-green)](LICENSE)
[![python: 3.11](https://img.shields.io/badge/python-3.11-44cc12?style=flat-square&logo=python)](https://www.python.org/downloads/release/python-3110/)
[![C++](https://img.shields.io/badge/C++-17-44cc12?style=flat-square&logo=c%2B%2B)](https://isocpp.org/)
[![docs: ](https://readthedocs.org/projects/ebonite/badge/?style=flat-square)](https://Globalizer.readthedocs.io/ru/latest/)
[![build:](https://github.com/UNN-ITMM-Software/iOpt/actions/workflows/python-app.yml/badge.svg)](https://github.com/OptimLLab/Globalizer/actions)



# Globalizer — программная система для параллельного поиска глобально-оптимальных решений
## Краткое описание
Программная система `Globalizer` предназначена для решения трудоемких задач многомерной многоэкстремальной глобальной оптимизации.

Методы глобальной оптимизации, реализованные в составе программной системы Globalizer, предназначены для выбора значений параметров математических моделей сложных объектов и процессов. Система Globalizer позволяет проводить точную настройку параметров моделей и методов, используемых в различных прикладных областях. Характерными примерами таких задач являются задачи настройки гиперпараметров методов ИИ и МО.

Параллельная версия алгоритма глобального поиска реализована с использованием `OpenMP`, `MPI`, `CUDA.

# **Основные возможности**
- Автоматический выбор значений параметров математических моделей и методов ИИ и МО, используемых в промышленности.
- Интеллектуальное управление процессом выбора оптимальных параметров для промышленных задач.
- Интеграция с внешними библиотеками или фреймворками искусственного интеллекта и машинного обучения, а также предметными моделями.
- Визуализация процесса выбора оптимальных параметров.


# **Установка и настройка**

## Сборка Globalizer со встроенными задачами:

```
git clone https://github.com/OptimLLab/Globalizer.git
cd Globalizer/
mkdir build
cd build
cmake -DGLOBALIZER_USE_MPI=OFF -DGLOBALIZER_USE_MP=OFF -DGLOBALIZER_USE_CUDA=OFF -DGLOBALIZER_BUILD_TESTS=OFF ..
cmake --build ../build --config Release
../_bin/GlobalizerSimpleMain.exe
```

## Сборка Globalizer для Microsoft Visual Studio 2022 (*У Вас должена быть установлена и выбрана средой по умолчанию Microsoft Visual Studio 2022):
```
git clone https://github.com/OptimLLab/Globalizer.git
cd Globalizer/gen
./vs-17_64-no_mp-No_MPI.bat
*Выберите запускаемым проектом GlobalizerSimpleMain соберите решение и запустите приложение*
```

## Сборка Globalizer с задачами из репозиториев https://github.com/OptimLLab/Globalizer_Benchmarks и https://github.com/OptimLLab/GCGen:

```
git clone https://github.com/OptimLLab/Globalizer.git
cd Globalizer/
mkdir build
cd build
conda init
conda create -p /Globalizer/Globalizer_env python=3.11
conda activate /Globalizer/Globalizer_env
conda install numpy
conda install pytorch::pytorch 
conda install conda-forge::pytorch-lightning
conda install lightning
conda install scikit-learn
git submodule init
git submodule update
cmake -DGLOBALIZER_USE_MPI=OFF -DGLOBALIZER_USE_MP=OFF -DGLOBALIZER_USE_CUDA=OFF -DGLOBALIZER_BUILD_TESTS=OFF -DGLOBALIZER_BUILD_PROBLEMS=ON -DBUILD_ALL_TASK=ON ..
cmake --build ../build --config Release
cd ../_bin/
GlobalizerSampleMain.exe -libPath rastrigin.dll -N 4
```


# **Начать работать**

Использование фреймворка Globalizer для минимизации функции Стронгина с тремя ограничениями.

```\examples\SimpleMain.cpp

#ifndef _CRT_SECURE_NO_WARNINGS
#define _CRT_SECURE_NO_WARNINGS
#endif

#include "Solver.h"
#include "GlobalizerProblem.h"

#include <algorithm>
#include <fstream>
#include <vector>
#include <string>
#include <cmath>
#include <iostream>

enum ProblemName { RASTRIGIN, STRONGINC3_LAMBDA_EXPRESSION, STRONGINC3_FUNCTION_POINTER};

double StronginC3Functionals(const double* y, int fNumber)
{
  double res = 0.0;
  double x1 = y[0], x2 = y[1];
  switch (fNumber)
  {
  case 0: // constraint 1
    res = 0.01 * ((x1 - 2.2) * (x1 - 2.2) + (x2 - 1.2) * (x2 - 1.2) - 2.25);
    break;
  case 1: // constraint 2
    res = 100.0 * (1.0 - ((x1 - 2.0) / 1.2) * ((x1 - 2.0) / 1.2) -
      (x2 / 2.0) * (x2 / 2.0));
    break;
  case 2: // constraint 3
    res = 10.0 * (x2 - 1.5 - 1.5 * sin(6.283 * (x1 - 1.75)));
    break;
  case 3: // criterion
  {
    double t1 = pow(0.5 * x1 - 0.5, 4.0);
    double t2 = pow(x2 - 1.0, 4.0);
    res = 1.5 * x1 * x1 * exp(1.0 - x1 * x1 - 20.25 * (x1 - x2) * (x1 - x2));
    res = res + t1 * t2 * exp(2.0 - t1 - t2);
    res = -res;
  }
  break;
  }

  return res;
}

// ------------------------------------------------------------------------------------------------
int main(int argc, char* argv[])
{
  parameters.Init(argc, argv, true);

  // Инициализация системы вывода и печати ошибок
  OutputMessage::Init(true, parameters.logFileNamePrefix, parameters.GetProcNum(),
    parameters.GetProcRank());


  parameters.Dimension = 2;
  ProblemName problemName = STRONGINC3_FUNCTION_POINTER;
  IProblem* problem = nullptr;

  if (problemName == RASTRIGIN)
  {
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      std::vector<double>(parameters.Dimension, -2.2), // нижняя граница
      std::vector<double>(parameters.Dimension, 1.8), //  верхняя граница
      std::vector<std::function<double(const double*)>>(1, [](const double* y)
        {
          double pi_ = 3.14159265358979323846;
          double sum = 0.;
          for (int j = 0; j < parameters.Dimension; j++)
            sum += y[j] * y[j] - 10. * cos(2.0 * pi_ * y[j]) + 10.0;
          return sum;
        }), // критерий
      true, // определен ли оптимум
      0, // значение глобального оптимума
      std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума
    );
  }
  else if (problemName == STRONGINC3_LAMBDA_EXPRESSION)
  {
    parameters.r = 4;
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      {0.0, -1.0}, // нижняя граница
      {4.0, 3.0}, // верхняя граница
      std::vector<std::function<double(const double*)>>({ 
        [](const double* y) { return 0.01 * ((y[0] - 2.2) * (y[0] - 2.2) + (y[1] - 1.2) * (y[1] - 1.2) - 2.25); }, // ограничение 0
        [](const double* y) { return 100.0 * (1.0 - ((y[0] - 2.0) / 1.2) * ((y[0] - 2.0) / 1.2) - (y[1] / 2.0) * (y[1] / 2.0)); }, // ограничение 1
        [](const double* y) { return 10.0 * (y[1] - 1.5 - 1.5 * sin(6.283 * (y[0] - 1.75))); }, // ограничение 2
        [](const double* y) 
        { 
          double t1 = pow(0.5 * y[0] - 0.5, 4.0);
          double t2 = pow(y[1] - 1.0, 4.0);
          return -((1.5 * y[0] * y[0] * exp(1.0 - y[0] * y[0] - 20.25 * (y[0] - y[1]) * (y[0] - y[1]))) + t1 * t2 * exp(2.0 - t1 - t2));
        } // ограничение 2
        }), 
      true, // определен ли оптимум
      0, // значение глобального оптимума
      std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума
    );
  }
  else if (problemName == STRONGINC3_FUNCTION_POINTER)
  {
    parameters.r = 4;
    problem = new ProblemFromFunctionPointers(parameters.Dimension, // размерность задачи
      { 0.0, -1.0 }, // нижняя граница
      { 4.0, 3.0 }, // верхняя граница
      StronginC3Functionals, // задача
      4, // количество функций (3 ограничения + 1 критерий)
      true, // определен ли оптимум
      0, // значение глобального оптимума
      std::vector<double>(parameters.Dimension, 0).data() // координаты глобального минимума
    );
  }

  problem->Initialize();

  // Решатель
  Solver solver(problem);

  // Решаем задачу
  if (solver.Solve() != SYSTEM_OK)
    throw EXCEPTION("Error: solver.Solve crash!!!");

  return 0;
}
```

# **Примеры**


# **Структура проекта**

Последняя стабильная версия Globalizer находится в [ветке main](https://github.com/OptimLLab/Globalizer/tree/main). 

Репозиторий включает в себя следующие директории::
- Пакет [globalizer](https://github.com/OptimLLab/Globalizer/tree/main/globalizer) содержит ядро фреймворка  в виде  классов на языке C++.
- Пакет [examples](https://github.com/OptimLLab/Globalizer/tree/main/examples) содержит примеры применения фреймворка для модельных задач.
- Модульные тесты размещены в каталоге [test](https://github.com/OptimLLab/Globalizer/tree/main/tests).
- Исходные файлы документации находятся в каталоге [docs](https://github.com/OptimLLab/Globalizer/tree/main/docs).

# **Документация**

Подробное описание API фреймворка Globalizer доступно в разделе [Read the Docs](https://Globalizer.readthedocs.io/ru/latest/).

# **Поддержка**

Исследование проводится при поддержке :