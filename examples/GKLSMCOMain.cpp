#include "Globalizer.h"
#include "MCOSolver.h"

#include <algorithm>
#include <array>
#include <charconv>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <limits>
#include <stdexcept>
#include <string>
#include <vector>

namespace
{
int PositiveInteger(const std::string& value)
{
  int result = 0;
  auto parsed = std::from_chars(value.data(), value.data() + value.size(), result);
  if (parsed.ec != std::errc() || parsed.ptr != value.data() + value.size() || result < 1)
    throw std::invalid_argument("Expected a positive integer");
  return result;
}

class CountedProblem : public GlobalizerBenchmarksProblem
{
public:
  using GlobalizerBenchmarksProblem::GlobalizerBenchmarksProblem;
  std::array<int, 4> calls{};
  double CalculateFunctionals(const double* x, int index) override
  {
    ++calls.at(index);
    return GlobalizerBenchmarksProblem::CalculateFunctionals(x, index);
  }
};

void Grid(IGlobalOptimizationProblem& problem, int size, bool map, std::ostream& out)
{
  if (size < 2)
    throw std::invalid_argument("Grid side must be at least two");
  std::vector<std::array<double, 4>> feasible;
  std::vector<double> x(2);
  std::vector<std::string> discrete;
  if (map)
    out << "x,y,g1,g2\n";
  for (int i = 0; i < size; ++i)
    for (int j = 0; j < size; ++j)
    {
      x[0] = -1.0 + 2.0 * i / (size - 1);
      x[1] = -1.0 + 2.0 * j / (size - 1);
      const double g1 = problem.CalculateFunctionals(x, discrete, 0);
      if (map)
        out << x[0] << ',' << x[1] << ',' << g1 << ',' << problem.CalculateFunctionals(x, discrete, 1) << '\n';
      else if (g1 <= 0 && problem.CalculateFunctionals(x, discrete, 1) <= 0)
        feasible.push_back({ problem.CalculateFunctionals(x, discrete, 2),
                             problem.CalculateFunctionals(x, discrete, 3), x[0], x[1] });
    }
  if (map)
    return;
  std::sort(feasible.begin(), feasible.end());
  double best = std::numeric_limits<double>::infinity();
  out << "f1,f2,x,y\n";
  for (const auto& p : feasible)
    if (p[1] < best)
    {
      out << p[0] << ',' << p[1] << ',' << p[2] << ',' << p[3] << '\n';
      best = p[1];
    }
  std::cout << "feasible=" << feasible.size() << '\n';
}

void Solve(IGlobalOptimizationProblem& problem, int budget, double rLocal, std::ostream& out)
{
  if (budget < 50 || !std::isfinite(rLocal) || (rLocal != 0 && rLocal <= 1))
    throw std::invalid_argument("Expected budget >= 50 and rLocal = 0 or > 1");
  parameters.Dimension = 2;
  parameters.Epsilon = 0.02;
  parameters.r = 5.6;
  parameters.m = 12;
  parameters.MaxNumOfPoints = budget;
  parameters.IsPrintResultToConsole = false;
  MCOParameters mco;
  mco.numberOfLambdas = 50;
  mco.rLocal = rLocal;
  CountedProblem wrapped(&problem);
  MCOSolver solver(&wrapped, mco);
  if (solver.Solve() != SYSTEM_OK)
    throw std::runtime_error("MCOSolver failed");
  out << "lambda,x,y,f1,f2\n";
  for (const auto& p : solver.GetTrials())
  {
    out << p.lambdaIndex << ',' << p.y[0] << ',' << p.y[1] << ',';
    if (p.criteria.empty())
      out << ',';
    else
      out << p.criteria[0] << ',' << p.criteria[1];
    out << '\n';
  }
  std::cout << "calls=" << wrapped.calls[0] << ',' << wrapped.calls[1] << ','
            << wrapped.calls[2] << ',' << wrapped.calls[3] << '\n';
}
}

int main(int argc, char** argv)
{
  if (argc < 6)
  {
    std::cerr << "Usage: GKLSMCOMain library instance output.csv grid|map side | solve budget rLocal\n";
    return 1;
  }
  try
  {
    GlobalizerInitialization(1, argv);
    GlobalOptimizationProblemManager manager;
    if (manager.LoadProblemLibrary(argv[1]) != 0)
      throw std::runtime_error("Cannot load benchmark library");
    auto* problem = manager.GetProblem();
    if (problem->SetParameter("problemIndex", std::string(argv[2])) != 0 ||
        problem->Initialize() != 0 || problem->GetDimension() != 2 ||
        problem->GetNumberOfConstraints() != 2 || problem->GetNumberOfCriterions() != 2)
      throw std::runtime_error("Expected a 2D GKLS MCO instance in [1,100]");
    std::ofstream out;
    out.exceptions(std::ios::badbit | std::ios::failbit);
    out.open(argv[3]);
    out << std::setprecision(17);
    const int limit = PositiveInteger(argv[5]);
    const std::string mode = argv[4];
    if ((mode == "grid" || mode == "map") && argc == 6)
      Grid(*problem, limit, mode == "map", out);
    else if (mode == "solve" && argc == 7)
    {
      const std::string value = argv[6];
      size_t used = 0;
      const double rLocal = std::stod(value, &used);
      if (used != value.size())
        throw std::invalid_argument("Invalid rLocal");
      Solve(*problem, limit, rLocal, out);
    }
    else
      throw std::invalid_argument("Invalid mode or number of arguments");
    out.close();
  }
  catch (const std::exception& e)
  {
    std::cerr << e.what() << '\n';
    return 1;
  }
  return 0;
}
