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
    Real sigma1 = getSigma1();
    Real sigma2 = getSigma2();
    Real dt = this->option_->getExpiry() / this->N_tau_;
    auto finite_diff_bs = [&dt](Real S, Real dS, Real sigma,
                             Real r, Real up, Real cen, Real down){
        return cen + dt *(
            ((sigma * sigma * S * S / 2) * (up - 2 * cen + down) / dS) +
            (r * S * (up - down) / (2 * dS)) -
            r * cen
        );
    };

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

    for (uint32_t t = 0; t < this->N_tau_ - 2; ++t)
    {
        for (uint32_t i = 1; i < N1_ - 1; ++i)
        {
            auto S1 = i * dS1_;
            auto up = U(t, (i + 1), 0);
            auto cen = U(t, i, 0);
            auto down = U(t, (i - 1), 0);
            U(t + 1, i, 0) = finite_diff_bs(S1, dS1_, sigma1, r, up, cen, down);
        }
    }

    for (uint32_t t = 0; t < this->N_tau_ - 2; ++t)
    {
        for (uint32_t j = 1; j < N2_ - 1; ++j)
        {
            auto S2 = j * dS2_;
            auto up = U(t, 0, (j + 1));
            auto cen = U(t, 0, j);
            auto down = U(t, 0, (j - 1));
            U(t + 1, 0, j) = finite_diff_bs(S2, dS2_, sigma2, r, up, cen, down);
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

template class HJBSolver<double>;
