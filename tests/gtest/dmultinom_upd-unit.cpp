#include "gtest/gtest.h"
#include "../../inst/include/AuxillaryFuns.hpp"

namespace {
// inspired/stolen from https://noaa-fims.github.io/collaborative_workflow/testing.html#c-unit-testing-and-benchmarking
// TestSuiteName: dmultinom_updTest; TestName: DoubleInput
// Test dmultinom_upd with double input values

TEST(dmultinom_updTest, DoubleInput) {
  /* R code to reproduce this
   * set.seed(123)
   * n = 10
   * prob = rlnorm(n, log(0.4), 0.2)
   * prob = prob / sum(prob)
   * x = rmultinom(1, size = n, prob = prob)
   * stats::dmultinom(x, prob = prob,log = T)
   *
   */
  vector<double> input_prob = {0.087, 0.093, 0.132, 0.098, 0.099, 0.136, 0.106, 0.075, 0.084, 0.089};
  vector<double> x = {2,    1,    1 ,   3    ,1  ,  1   , 0   , 0  ,  0,     1};

  EXPECT_NEAR( dmultinom_upd(x, input_prob, 1) , -10.3614, 0.0001 );

}
