#include <gtest/gtest.h>
#include <functional>
#include "Option.hpp"
#include "PDESolver.hpp"


TEST(PDESolverTest, BasicPDESolverAssertion)
{
    std::function<double(double, double)> payoff_fn = [](double S1, double S2){
        return std::max(S1 - S2, 0.0);
    };
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.1, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.1, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    Option<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    size_t N_tau = 25;
    size_t N_12 = 91;

    PDESolver<double> pde_solver(N_12, N_12, N_tau, opt);
    auto U = pde_solver.getU();

    EXPECT_EQ(U(0, 0, 0), 0.0);
}