#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <vector>
#include <functional>
#include <concepts>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <Option.hpp>

constexpr uint16_t T_DIM = 3;

// ! Numerical HJB PDE Solver Class
template<std::floating_point Real>
class PDESolver {
    public:

        PDESolver(const uint32_t N1, const uint32_t N2, const uint32_t N_tau,
                  Option<Real>& option, std::pair<Real, Real>& S_max, bool is_sup);

        void solve();

        Eigen::Tensor<Real, T_DIM>& getU();
    
    protected:
        Eigen::Tensor<Real, T_DIM> initBoundaryConditions();

        Real getSigma1();

        Real getSigma2();

        const uint32_t N1_, N2_, N_tau_;
        Eigen::Tensor<Real, T_DIM> U_;
        Option<Real> option_;
        const double dS1_, dS2_;
        const std::pair<Real, Real> S_max_;
        bool is_sup_;
};

template<std::floating_point Real>
class FixedStencil : public PDESolver<Real> {
    public:
        // Constructor
        FixedStencil();

        void solve() override;

    private:

        double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2);
        void solve_sparse_matrix(std::vector<std::vector<double>>& U);
        void policy_iteration(std::vector<std::vector<double>>& U);
        void apply_boundary_conditions(std::vector<std::vector<double>>& U);
};


#endif //HH__PDESOLVERSHPP_HH