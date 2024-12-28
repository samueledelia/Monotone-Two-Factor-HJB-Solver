#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <gtest/gtest.h>
#include "BSSolver.hpp"
#include "HJBSolver.hpp"
#include "Option.hpp"
#include "TestUtils.hpp"

const std::string test_data_dir = TEST_DATA_DIR;

constexpr double TOL = 1e-2;

TEST(PDESolverTest, BasicPDESolverAssertion)
{
    std::function payoff_fn = [](double S1, double S2){
        double K = 100;
        return std::max(std::max(S1, S2) - K, 0.0);
    };
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.1, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.1, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    auto opt = std::make_unique<TwoAssetMinMaxOption<double>>(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    const uint32_t N_tau = 25;
    const uint32_t N_12 = 91;
    std::pair<double, double> S_max = std::make_pair(400, 400); 

    HJBSolver pde_solver(N_12, N_12, N_tau, std::move(opt), S_max, true);
    auto U = pde_solver.getU();

    Eigen::array<int, 2> shuffling({1, 0});
    Eigen::Tensor<double, 2> expected_sol_BC1 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/left_bound_hjb.txt");
    Eigen::Tensor<double, 2> expected_sol_BC2 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/right_bound_hjb.txt");
    Eigen::Tensor<double, 2> expected_sol_BC3 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/upper_bound_hjb.txt");
    Eigen::Tensor<double, 2> expected_sol_BC4 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/lower_bound_hjb.txt");
    Eigen::Tensor<double, 2> expected_sol_BC5 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/front_bound_hjb.txt");
    Eigen::Tensor<double, 2> expected_sol_BC6 = readMatrixFromFile<Eigen::Tensor<double, 2>>(test_data_dir + "/behind_bound_hjb.txt");
    Eigen::Tensor<double, 2> left_bound = U.chip(0, 2).shuffle(shuffling);
    Eigen::Tensor<double, 2> right_bound = U.chip(U.dimension(2) - 1, 2).shuffle(shuffling);
    Eigen::Tensor<double, 2> upper_bound = U.chip(U.dimension(1) - 1, 1).shuffle(shuffling);
    Eigen::Tensor<double, 2> lower_bound = U.chip(0, 1).shuffle(shuffling);
    Eigen::Tensor<double, 2> front_bound = U.chip(0, 0).shuffle(shuffling);
    Eigen::Tensor<double, 2> behind_bound = U.chip(U.dimension(0) - 1, 0).shuffle(shuffling);
    // Check boundary conditions
    EXPECT_EQ(U(N_tau - 1, 0, 0), 0.0);
    EXPECT_EQ(U(N_tau - 1, 23, 22), 1.0989010989011092);
    EXPECT_EQ(U(N_tau - 1, 23, 0), 1.0989010989011092);
        //<< std::fixed << std::setprecision(2) << left_bound;
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC1, left_bound, TOL));
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC2, right_bound, TOL));
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC3, upper_bound, TOL));
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC4, lower_bound, TOL));
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC5, front_bound, TOL));
    EXPECT_TRUE(TensorAreApproxEqual(expected_sol_BC6, behind_bound, TOL));
}

TEST(BSSolverTest, BasicBSSolverAssertion)
{
    std::function payoff_fn = [](double S1){
        double K = 100;
        return std::max(S1 - K, 0.0);
    };
    double r = 0.05;
    const auto sigma = 0.3;
    double expiry = 1.0;

    auto opt = std::make_unique<OneAssetOption<double>>(payoff_fn, r, sigma, expiry);

    const uint32_t N_tau = 25;
    const uint32_t N_12 = 91;
    double S_max = 200;

    BSSolver pde_solver(N_12, N_tau, std::move(opt), S_max);
    pde_solver.solve();
    auto U = pde_solver.getV();

    auto expected_sol = readMatrixFromFile<Eigen::MatrixXd>(test_data_dir + "/bs_pde_sol_1.txt");

    EXPECT_TRUE(expected_sol.isApprox(U, TOL))
        << "Calculated matrix does not match expected matrix";
}