#include <gtest/gtest.h>
#include <functional>
#include "Option.hpp"
#include "PDESolver.hpp"


TEST(PDESolverTest, BasicPDESolverAssertion)
{
    std::function<double(double, double)> payoff_fn = [](double S1, double S2){
        double K = 5;
        return std::max(std::max(S1 - S2, 0.0) - K, 0.0);
    };
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.1, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.1, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    Option<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    const uint32_t N_tau = 25;
    const uint32_t N_12 = 91;
    std::pair<double, double> S_max = std::make_pair(400, 400); 

    PDESolver<double> pde_solver(N_12, N_12, N_tau, opt, S_max, true);
    auto U = pde_solver.getU();

    // Check Boundary
    EXPECT_EQ(U(0, 0, 0), 0.0);
    EXPECT_EQ(U(0, 41, 39), 3.7912087912087884); // Coincide with Payoff = Max[(41-39)*dS -K, 0]
    EXPECT_EQ(U(5, 2, 0), 3.97990); // B&S Closed formula For t = 5 * (expiry/N_tau), S1 = 2 * dS

}