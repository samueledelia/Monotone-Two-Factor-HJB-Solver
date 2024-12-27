#include "HJBSolver.hpp"
#include <BSSolver.hpp>

template <std::floating_point Real>
HJBSolver<Real>::HJBSolver(const uint32_t N1, const uint32_t N2, const uint32_t N_tau,
    std::unique_ptr<TwoAssetMinMaxOption<Real>> option,
    std::pair<Real, Real>& S_max, bool is_sup)
    : PDESolver<Real, HJB_T_DIM>(std::move(option), N_tau), N1_(N1), N2_(N2),
    dS1_(S_max.first / N1), dS2_(S_max.second / N2), S_max_(S_max), is_sup_(is_sup)
{   
    S1_ = (Vector::LinSpaced(N1, 0, S_max.first) * dS1_).array();
    S2_ = (Vector::LinSpaced(N2, 0, S_max.second) * dS2_).array();
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
Real HJBSolver<Real>::getU_element(size_t i, size_t j) const {
    if (i < static_cast<size_t>(this->U_.dimension(1)) && j < static_cast<size_t>(this->U_.dimension(2))) {
        return this->U_(0, i, j);
    } else {
        throw std::out_of_range("Index out of range in getU_element.");
    }
}


template <std::floating_point Real>
Real HJBSolver<Real>::getS1(size_t i) const {
    if (static_cast<typename Eigen::EigenBase<Eigen::Matrix<Real, -1, 1>>::Index>(i) < S1_.size()) {
        return S1_[i];
    } else {
        throw std::out_of_range("Index out of range in getS1.");
    }
}

template <std::floating_point Real>
Real HJBSolver<Real>::getS2(size_t i) const {
    if (static_cast<typename Eigen::EigenBase<Eigen::Matrix<Real, -1, 1>>::Index>(i) < S2_.size()) {
        return S2_[i];
    } else {
        throw std::out_of_range("Index out of range in getS2.");
    }
}


// Helper to compute the cross-derivative term Γ
template <std::floating_point Real>
Real HJBSolver<Real>::computeGamma(Real rho, uint32_t i, uint32_t j) const
{
    // Retrieve the values of U at the necessary indices
    Real Uni_j = this->getU_element(i, j);         // U(i,j)
    Real Uni1_j1 = this->getU_element(i+1, j+1);   // U(i+1, j+1)
    Real Uni1_j = this->getU_element(i+1, j);      // U(i+1, j)
    Real Uni_j1 = this->getU_element(i, j+1);      // U(i, j+1)
    Real Uni_m1_j_m1 = this->getU_element(i-1, j-1); // U(i-1, j-1)
    //Real Uni_m1_j1 = this->getU_element(i-1, j+1); // U(i-1, j+1)
    Real Uni_m1_j = this->getU_element(i-1, j);    // U(i-1, j)
    Real Uni_j_m1 = this->getU_element(i, j-1);    // U(i, j-1)
    Real Uni1_j_m1 = this->getU_element(i+1, j-1); // U(i+1, j-1)
    

    // Compute Δ+(S1)i, Δ−(S1)i, Δ+(S2)j, Δ−(S2)j
    Real delta_plus_S1 = this->getS1(i+1) - this->getS1(i); // Δ+(S1)i
    Real delta_minus_S1 = this->getS1(i) - this->getS1(i-1); // Δ−(S1)i
    Real delta_plus_S2 = this->getS2(j+1) - this->getS2(j); // Δ+(S2)j
    Real delta_minus_S2 = this->getS2(j) - this->getS2(j-1); // Δ−(S2)j

    // Compute the cross-derivative term based on the sign of rho
    if (rho >= 0) {
        // Formula (4.4) for ρ >= 0
        Real term1 = (2 * Uni_j + Uni1_j1 + Uni_m1_j_m1) / (delta_plus_S1 * delta_plus_S2 + delta_minus_S1 * delta_minus_S2);
        Real term2 = (Uni1_j + Uni_m1_j + Uni_j1 + Uni_j_m1) / (delta_plus_S1 * delta_plus_S2 + delta_minus_S1 * delta_minus_S2);
        return term1 - term2;
    } else {
        // Formula (4.5) for ρ < 0
        Real term1 = -(2 * Uni_j + Uni1_j_m1 + Uni_m1_j_m1) / (delta_plus_S1 * delta_minus_S2 + delta_minus_S1 * delta_plus_S2);
        Real term2 = (Uni1_j + Uni_m1_j + Uni_j1 + Uni_j_m1) / (delta_plus_S1 * delta_minus_S2 + delta_minus_S1 * delta_plus_S2);
        return term1 + term2;
    }
}

/*
template <std::floating_point Real>
std::vector HJBSolver<Real>::determineOptimalControl(
    const Eigen::Tensor<Real, HJB_T_DIM>& W, 
    uint32_t i, uint32_t j) 
{
    // Step 1: Determine the optimal ρ
    Real rho_min = this->option_->getRhos().first;
    Real rho_max = this->option_->getRhos().second;

    Real gamma_min = this->computeGamma(rho_min, i, j);
    Real gamma_max = this->computeGamma(rho_max, i, j);

    Real rho_opt = (rho_max * gamma_max >= rho_min * gamma_min) ? rho_max : rho_min;

    // Step 2: Compute the positive coefficient sets for central and forward differencing
    auto [sigma_central, sigma_forward] = this->computePositiveCoefficientSets(i, j);

    // Step 3: Initialize variables for tracking optimal values
    Real F_max = -std::numeric_limits<Real>::infinity();
    Control optimal_control{0.0, 0.0, rho_opt}; // (sigma1, sigma2, rho)

    // Step 4: Iterate over differencing methods
    for (auto& differencing : {sigma_central, sigma_forward}) 
    {
        for (auto& interval : differencing) 
        {
            // Solve the quadratic problem on each interval
            auto [sigma1_opt, sigma2_opt] = this->maximizeObjectiveFunction(interval, W, rho_opt);

            // Evaluate the objective function
            Real F_value = this->evaluateObjectiveFunction(sigma1_opt, sigma2_opt, rho_opt, W);

            // Check if this is the new maximum
            if (F_value > F_max) 
            {
                F_max = F_value;
                optimal_control = {sigma1_opt, sigma2_opt, rho_opt};
            }
        }
    }

    return optimal_control;
}

*/
template <std::floating_point Real>
void HJBSolver<Real>::solve(){

}
template class HJBSolver<double>;