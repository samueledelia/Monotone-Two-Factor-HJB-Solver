#include <gtest/gtest.h>
#include <functional>
#include "Option.hpp"
#include "HJBSolver.hpp"
#include "BSSolver.hpp"
#include <iomanip>


TEST(PDESolverTest, BasicPDESolverAssertion)
{
    std::function<double(double, double)> payoff_fn = [](double S1, double S2){
        double K = 100;
        return std::max(std::max(S1, S2) - K, 0.0);
    };
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.1, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.1, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    Option<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    const uint32_t N_tau = 252;
    const uint32_t N_12 = 91;
    std::pair<double, double> S_max = std::make_pair(400, 400); 

    HJBSolver<double> pde_solver(N_12, N_12, N_tau, opt, S_max, true);
    auto U = pde_solver.getU();

    // Check boundary conditions
    EXPECT_EQ(U(N_tau - 1, 0, 0), 0.0);
    EXPECT_EQ(U(N_tau - 1, 23, 22), 1.0989010989011092);
    EXPECT_EQ(U(N_tau - 1, 23, 0), 1.0989010989011092);
        //<< std::fixed << std::setprecision(2) << U.chip(0, 2).chip(23, 1);
}

TEST(BSSolverTest, BasicBSSolverAssertion)
{
    std::function<double(double, double)> payoff_fn = [](double S1,[[maybe_unused]] double S2){
        double K = 100;
        return std::max(S1 - K, 0.0);
    };
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.3, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.3, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    Option<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    const uint32_t N_tau = 25;
    const uint32_t N_12 = 91;
    double S_max = 200;

    BSSolver<double> pde_solver(N_12, N_tau, opt, S_max);
    pde_solver.solve();
    auto U = pde_solver.getV();

    // Check boundary conditions
    EXPECT_EQ(U( 46, N_tau - 1), 2.2222222222222143);
        //<< std::fixed << std::setprecision(2) << U;
}