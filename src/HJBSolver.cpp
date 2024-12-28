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

    //U(0, 0, 0) = std::exp(-r) * this->option_->evaluate(0, 0);

    // Boundary condition for S2 = 0
    auto payoff = [this](Real S1) { return this->option_->evaluate(S1, 0); };
    auto bs_option = std::make_unique<OneAssetOption<Real>>(payoff, r, this->option_->getSigmas_1().second, this->option_->getExpiry()); //GETSIGMA MIN OR MAX?
    BSSolver<Real> bs_solver(N1_, this->N_tau_, std::move(bs_option), S_max_.first);

    bs_solver.solve();
    Eigen::Tensor V = bs_solver.getU();

    for (uint32_t t = 0; t < this->N_tau_ - 1; ++t)
    {
        for (uint32_t i = 0; i < N1_ ; ++i)
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

    for (uint32_t t = 0; t < this->N_tau_ - 1; ++t)
    {
        for (uint32_t j = 0; j < N2_ ; ++j)
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

    for (int32_t t = this->N_tau_ - 1; t >= 0; t--)
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

/*
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
        // for ρ >= 0
        Real term1 = (2 * Uni_j + Uni1_j1 + Uni_m1_j_m1) / (delta_plus_S1 * delta_plus_S2 + delta_minus_S1 * delta_minus_S2);
        Real term2 = (Uni1_j + Uni_m1_j + Uni_j1 + Uni_j_m1) / (delta_plus_S1 * delta_plus_S2 + delta_minus_S1 * delta_minus_S2);
        return term1 - term2;
    } else {
        // for ρ < 0
        Real term1 = -(2 * Uni_j + Uni1_j_m1 + Uni_m1_j_m1) / (delta_plus_S1 * delta_minus_S2 + delta_minus_S1 * delta_plus_S2);
        Real term2 = (Uni1_j + Uni_m1_j + Uni_j1 + Uni_j_m1) / (delta_plus_S1 * delta_minus_S2 + delta_minus_S1 * delta_plus_S2);
        return term1 + term2;
    }
}

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




template <std::floating_point Real>
void HJBSolver<Real>::buildLinearSystem(
    Eigen::SparseMatrix<Real>& A, 
    Eigen::Vector<Real, Eigen::Dynamic>& C, 
    const Eigen::Vector<Real, Eigen::Dynamic>& U_guess, 
    uint32_t t) const 
{
    // Dimensions
    const size_t nS1 = this->S1_.size();
    const size_t nS2 = this->S2_.size();
    const size_t nVars = nS1 * nS2;

    // Resize matrix and vector
    A.resize(nVars, nVars);
    C.resize(nVars);

    // Constants and time step
    const Real dt = this->dt_;
    const Real sigma1 = this->sigma1_;
    const Real sigma2 = this->sigma2_;
    const Real rho = this->rho_;
    const Real r = this->r_;
    const Real q1 = this->q1_;
    const Real q2 = this->q2_;

    // Populate A and C
    std::vector<Eigen::Triplet<Real>> triplets; // For efficient sparse matrix assembly

    for (size_t j = 0; j < nS2; ++j) {
        for (size_t i = 0; i < nS1; ++i) {
            size_t index = i + j * nS1; // Single index notation for (i, j)
            Real L_value = 0.0;

            // Compute contributions from LQ_U
            if (i > 0 && i < nS1 - 1 && j > 0 && j < nS2 - 1) {
                // Internal nodes: apply central differencing
                L_value = LQ_U(this->U_, this->S1_, this->S2_, sigma1, sigma2, rho, r, q1, q2, i, j);
            } else {
                // Boundary conditions
                L_value = U_guess[index]; // Example: Dirichlet condition
            }

            // Construct A(Q) = [I - Δτ * L(Q)]
            triplets.emplace_back(index, index, 1.0 - dt * L_value);

            // Construct C(Q) = U^n
            C[index] = U_guess[index];
        }
    }

    // Finalize sparse matrix assembly
    A.setFromTriplets(triplets.begin(), triplets.end());
    A.makeCompressed();
}



template <std::floating_point Real>
void HJBSolver<Real>::solve() {
    // Constants
    const uint32_t max_iter = 100;  
    const Real tol = 1e-6;        

    // Initialize U within the solve method or constructor
    Eigen::Tensor<Real, HJB_T_DIM> U = this->initBoundaryConditions();

    // U_0 = U_payoff
    U_old = U.chip(this->N_tau_ - 1, 0);

    // Time-stepping loop
    for (uint32_t t = N_tau_ -2; t >= 0 ; --t) {

        // Flatten tensors into vectors
        auto flatten = [](const Eigen::Tensor<Real, HJB_T_DIM>& tensor) -> Eigen::Vector<Real, Eigen::Dynamic> {
            return Eigen::Map<const Eigen::Vector<Real, Eigen::Dynamic>>(tensor.data(), tensor.size());
        };
        auto unflatten = [](const Eigen::Vector<Real, Eigen::Dynamic>& vec, int N_tau, int N1, int N2) -> Eigen::Tensor<Real, HJB_T_DIM> {
            return Eigen::TensorMap<const Eigen::Tensor<Real, HJB_T_DIM>>(vec.data(), N_tau, N1, N2);
        };

        Eigen::Vector<Real, Eigen::Dynamic> U_vec = flatten(U_old);
        Eigen::Vector<Real, Eigen::Dynamic> U_guess = U_vec;

        // Policy iteration
        uint32_t iter = 0;
        bool converged = false;
        while (iter < max_iter && !converged) {
            // Build sparse matrix A(Q_k) and vector C(Q_k)
            Eigen::SparseMatrix<Real> A(U_guess.size(), U_guess.size());
            Eigen::Vector<Real, Eigen::Dynamic> C(U_guess.size());
            this->buildLinearSystem(A, C, U_guess, t); // Fill A and C based on the discrete equations

            // Solve A * U_vec = C using Bi-CGSTAB with preconditioning
            Eigen::BiCGSTAB<Eigen::SparseMatrix<Real>> solver;
            solver.setTolerance(tol);
            solver.setMaxIterations(1000);

            // Set preconditioner (e.g., Incomplete LU)
            Eigen::IncompleteLUT<Real> preconditioner;
            solver.compute(A);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Failed to factorize matrix A during preconditioning.");
            }

            U_vec = solver.solve(C);
            if (solver.info() != Eigen::Success) {
                throw std::runtime_error("Failed to solve linear system using Bi-CGSTAB.");
            }

            // Check convergence
            U_new = unflatten(U_vec, N_tau_, N1_, N2_);
            converged = (U_new - U_old).abs().maxCoeff() < tol;

            U_old = U_new;
            iter++;
        }

        // Check if policy iteration converged
        if (converged) {
            std::cout << "Time step " << t << ": converged in " << iter << " iterations." << std::endl;
        } else {
            std::cerr << "Time step " << t << ": did not converge within max iterations." << std::endl;
        }

        // Update solution for the next time step
        this->U_ = U_new;
    }
    
}
*/
template class HJBSolver<double>;