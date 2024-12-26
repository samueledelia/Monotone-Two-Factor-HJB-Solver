#include "BSSolver.hpp"

using SparseMatrix = Eigen::SparseMatrix<double>;
using Triplet = Eigen::Triplet<double>;

template<std::floating_point Real>
BSSolver<Real>::BSSolver(const uint32_t N1, const uint32_t N_tau,
                  std::unique_ptr<OneAssetOption<Real>> option, Real S_max)
    : PDESolver<Real, BS_T_DIM>(N_tau, std::move(option)), N_(N1), dS_(S_max / N1),
    dt_(1.0 / N_tau), S_max_(S_max)
{
    //this->U_ = BSSolver<Real>::initBoundaryConditions();
    S_ = (Eigen::VectorXd::LinSpaced(N_, 0, N_) * dS_).array();
    V_ = initBoundaryConditions_();
}

template <std::floating_point Real>
void BSSolver<Real>::solve()
{
    double sigma = getSigma();
    double r = this->option_->getDiscountRate();

    Vector I = Eigen::VectorXd::LinSpaced(N_, 0, N_);
    Eigen::VectorXd alpha = 0.25 * dt_ * ((sigma * sigma * I.array().square()).matrix() - r * I);
    Eigen::VectorXd beta = -0.5 * dt_ * ((sigma * sigma * I.array().square()).matrix() + Eigen::VectorXd::Constant(N_, r));
    Eigen::VectorXd gamma = 0.25 * dt_ * ((sigma * sigma * I.array().square()).matrix() + r * I);

    // Build sparse matrices ML and MR
    std::vector<Triplet> ML_triplets, MR_triplets;
    for (uint32_t i = 1; i < N_; ++i) {
        if (i > 1) {
            ML_triplets.emplace_back(i - 1, i - 2, -alpha(i));
            MR_triplets.emplace_back(i - 1, i - 2, alpha(i));
        }
        ML_triplets.emplace_back(i - 1, i - 1, 1.0 - beta(i));
        MR_triplets.emplace_back(i - 1, i - 1, 1.0 + beta(i));
        if (i < N_ - 1) {
            ML_triplets.emplace_back(i - 1, i, -gamma(i));
            MR_triplets.emplace_back(i - 1, i, gamma(i));
        }
    }

    SparseMatrix ML(N_ - 1, N_ - 1), MR(N_ - 1, N_ - 1);
    ML.setFromTriplets(ML_triplets.begin(), ML_triplets.end());
    MR.setFromTriplets(MR_triplets.begin(), MR_triplets.end());

    Eigen::SparseLU<SparseMatrix> solver;
    solver.compute(ML);
    if (solver.info() != Eigen::Success) {
        throw std::runtime_error("Matrix decomposition failed");
    }

    // Time-stepping loop
    for (int32_t t = static_cast<int32_t>(this->N_tau_) - 2; t >= 0; --t) {
        Eigen::VectorXd boundary_t = Eigen::VectorXd::Zero(N_ - 1);
        boundary_t(0) = alpha(1) * (V_(0, t) + V_(0, t + 1)) - alpha(0) * V_(0, t + 1);
        boundary_t(N_ - 2) = gamma(N_ - 1) * (V_(N_ - 1, t) + V_(N_ - 1, t + 1));

        Eigen::VectorXd b = MR * V_.block(1, t + 1, N_ - 1, 1) + boundary_t;
        Eigen::VectorXd x = solver.solve(b);
        if (solver.info() != Eigen::Success) {
            throw std::runtime_error("Solving failed");
        }
        V_.block(1, t, N_ - 1, 1) = x;
    }
}

template <std::floating_point Real>
Eigen::MatrixXd BSSolver<Real>::getV()
{
    return V_;
}

template <std::floating_point Real>
Eigen::Tensor<Real, BS_T_DIM> BSSolver<Real>::initBoundaryConditions()
{
    Eigen::Tensor<Real, BS_T_DIM> U(N_, this->N_tau_);
    U.setZero();
    double r = this->option_->getDiscountRate();
    double expiration = this->option_->getExpiry();

    // Initialize payoff and boundary conditions
    for (uint32_t i = 0; i < N_; ++i) {
        U(i, this->N_tau_ - 1) = this->option_->evaluate(S_(i));
    }
    for (uint32_t t = 0; t <= this->N_tau_; ++t) {
        double time = t * dt_;
        U(N_ - 1, t) = this->option_->evaluate(S_max_) * std::exp(-r * (expiration - time));
        U(0, t) = this->option_->evaluate(0.0) * std::exp(-r * (expiration - time));
    }

    return U;
}

template <std::floating_point Real>
Eigen::MatrixXd BSSolver<Real>::initBoundaryConditions_()
{
    Eigen::MatrixXd V = Eigen::MatrixXd::Zero(N_, this->N_tau_);
    double r = this->option_->getDiscountRate();
    double expiration = this->option_->getExpiry();

    // Initialize payoff and boundary conditions
    for (uint32_t i = 0; i < N_; ++i) {
        V(i, this->N_tau_ - 1) = this->option_->evaluate(S_(i));
    }
    for (uint32_t t = 0; t < this->N_tau_; ++t) {
        double time = t * dt_;
        V(N_ - 1, t) = this->option_->evaluate(S_max_) * std::exp(-r * (expiration - time));
        V(0, t) = this->option_->evaluate(0.0) * std::exp(-r * (expiration - time));
    }

    return V;
}

template <std::floating_point Real>
Real BSSolver<Real>::getSigma()
{
    return this->option_->getSigma();
}

template class BSSolver<double>;
