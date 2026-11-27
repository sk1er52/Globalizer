#ifdef GKLS_MCO_LIBRARY
#include <gtest/gtest.h>
#include "GlobalOptimizationProblemManager.h"
#include <array>

TEST(GKLSMultiCriteria, ReproducesLebedevFunctionValues)
{
  GlobalOptimizationProblemManager manager;
  ASSERT_EQ(0, manager.LoadProblemLibrary(GKLS_MCO_LIBRARY));
  auto* problem = manager.GetProblem();
  ASSERT_EQ(2, problem->GetNumberOfConstraints());
  ASSERT_EQ(2, problem->GetNumberOfCriterions());
  ASSERT_NE(0, problem->SetDimension(3));
  ASSERT_NE(0, problem->SetParameter("problemIndex", std::string("101")));

  // Значения g1, g2, f1, f2 в центре квадрата из реализации cbf4284.
  const std::array<int, 3> ids{1, 50, 100};
  const double expected[][4] = {
    {-2.0345374151789022, -1.3226670966896088, 0.93829319930198463, 0.10601220435387144},
    {-2.1453262905432515, -1.6630496141577389, 0.5452708189293024, 0.58037594567442141},
    {-0.042884688767597545, -0.018652320925943933, 0.62741569494035421, 0.93829319930198463}
  };
  std::vector<std::string> discrete;
  for (size_t i = 0; i < ids.size(); ++i)
  {
    ASSERT_EQ(0, problem->SetParameter("problemIndex", std::to_string(ids[i])));
    ASSERT_EQ(0, problem->Initialize());
    for (int j = 0; j < 4; ++j)
      EXPECT_NEAR(expected[i][j], problem->CalculateFunctionals({0., 0.}, discrete, j), 1e-11);
  }
}
#endif
