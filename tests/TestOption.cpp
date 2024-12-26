#include <gtest/gtest.h>
#include <functional>
#include "Option.hpp"

double payoff(double S1, double S2)
{
    return std::max(S1 - S2, 0.0);
}

TEST(OptionTest, BasicTwoAssetOptionAssertion)
{
    std::function<double(double, double)> payoff_fn = payoff;
    double r = 0.05;
    const auto sigmas_1 = 0.2;
    const auto sigmas_2 = 0.2;
    const auto rhos = 1.0;
    double expiry = 1.0;

    TwoAssetOption<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    EXPECT_EQ(opt.getExpiry(), expiry);
    EXPECT_EQ(opt.evaluate(110.0, 100.0), 10.0);
}

TEST(OptionTest, BasicTwoAssetMinMaxOptionAssertion)
{
    std::function<double(double, double)> payoff_fn = payoff;
    double r = 0.05;
    const auto sigmas_1 = std::make_pair<double, double>(0.1, 0.3);
    const auto sigmas_2 = std::make_pair<double, double>(0.1, 0.3);
    const auto rhos = std::make_pair<double, double>(-1.0, 1.0);
    double expiry = 1.0;

    TwoAssetMinMaxOption<double> opt(payoff_fn, r, sigmas_1, sigmas_2, rhos, expiry);

    EXPECT_EQ(opt.getExpiry(), expiry);
}