#include "HJBSolver.hpp"
#include <BSSolver.hpp>

template <std::floating_point Real>
HJBSolver<Real>::HJBSolver(const uint32_t N1, const uint32_t N2, const uint32_t N_tau,
    std::unique_ptr<TwoAssetMinMaxOption<Real>> option,
    std::pair<Real, Real>& S_max, bool is_sup)
    : PDESolver<Real, HJB_T_DIM>(std::move(option), N_tau), N1_(N1), N2_(N2),
    dS1_(S_max.first / N1), dS2_(S_max.second / N2), S_max_(S_max), is_sup_(is_sup)
{
    this->U_ = HJBSolver<Real>::initBoundaryConditions();
}

template <std::floating_point Real>
Eigen::Tensor<Real, HJB_T_DIM> HJBSolver<Real>::initBoundaryConditions()
{
    Eigen::Tensor<Real, HJB_T_DIM> U(this->N_tau_, N1_, N2_);
    U.setZero();
    Real r = this->option_->getDiscountRate();

    // Boundary condition at expiry
    for (uint32_t i = 0; i < N1_; ++i)
    {
        for (uint32_t j = 0; j < N2_; ++j)
        {
            auto S1 = i * dS1_;
            auto S2 = j * dS2_;
            U(this->N_tau_ - 1, i, j) = this->option_->evaluate(S1, S2);
        }
    }

    U(0, 0, 0) = std::exp(-r) * this->option_->evaluate(0, 0);

    // Boundary condition for S2 = 0
    auto payoff = [this](Real S1) { return this->option_->evaluate(S1, 0); };
    auto bs_option = std::make_unique<OneAssetOption<Real>>(payoff, r, this->option_->getSigmas_1().second, this->option_->getExpiry()); //GETSIGMA MIN OR MAX?
    BSSolver<Real> bs_solver(N1_, this->N_tau_, std::move(bs_option), S_max_.first);

    bs_solver.solve();
    Eigen::Tensor V = bs_solver.getU();

    for (uint32_t t = 0; t < this->N_tau_ - 2; ++t)
    {
        for (uint32_t i = 0; i < N1_ - 1; ++i)
        {
            U(t, i, 0) = V(i, t);
        }
    }

    // Boundary condition for S1 = 0
    auto payoff_S2 = [this](Real S2) { return this->option_->evaluate(0, S2); };
    auto bs_option_S2 = std::make_unique<OneAssetOption<Real>>(payoff_S2, r, this->option_->getSigmas_2().second, this->option_->getExpiry()); //GETSIGMA MIN OR MAX?
    BSSolver<Real> bs_solver_S2(N2_, this->N_tau_, std::move(bs_option_S2), S_max_.second);

    bs_solver_S2.solve();
    Eigen::Tensor V_S2 = bs_solver_S2.getU();

    for (uint32_t t = 0; t < this->N_tau_ - 2; ++t)
    {
        for (uint32_t j = 0; j < N2_ - 1; ++j)
        {
            U(t, 0, j) = V_S2(j, t);
        }
    }

    // Solve the ODEs and update upper(S1_max and S2_max) boundary values
    // Correct initial conditions for c0, c1, c2
    Real c0 = this->option_->evaluate(S_max_.first, S_max_.second);
    Real c1 = (this->option_->evaluate(S_max_.first, 0) - c0) / S_max_.first;
    Real c2 = (this->option_->evaluate(0, S_max_.second) - c0) / S_max_.second;
    Real dt = this->option_->getExpiry() / this->N_tau_;

    for (uint32_t t = this->N_tau_ - 2; t > 0; --t)
    {
        c0 *= std::exp(-r * dt);

        for (uint32_t i = 0; i < N1_; ++i)
        {
            auto S1 = i * dS1_;
            U(t, i, N2_ - 1) = c0 + c1 * S1 + c2 * S_max_.second;
        }
        for (uint32_t j = 0; j < N2_; ++j)
        {
            auto S2 = j * dS2_;
            U(t, N1_ - 1, j) = c0 + c1 * S_max_.first + c2 * S2;
        }
    }

    return U;
}

template <std::floating_point Real>
Real HJBSolver<Real>::getSigma1()
{
    return is_sup_ ? this->option_->getSigmas_1().second : this->option_->getSigmas_1().first;
}

template <std::floating_point Real>
Real HJBSolver<Real>::getSigma2()
{
    return is_sup_ ? this->option_->getSigmas_2().second : this->option_->getSigmas_2().first;
}

template <std::floating_point Real>
void HJBSolver<Real>::solve(){

}
template class HJBSolver<double>;