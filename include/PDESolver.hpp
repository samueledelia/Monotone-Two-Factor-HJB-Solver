#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <vector>
#include <functional>
#include <concepts>
#include <eigen3/Eigen/Dense>
#include <eigen3/unsupported/Eigen/CXX11/Tensor>
#include <Option.hpp>

// ! Numerical HJB PDE Solver Class
template<std::floating_point Real>
class PDESolver {
    public:

        PDESolver(size_t N1, size_t N2, size_t N_tau, Option<Real>& option);

        void solve();

        Eigen::Tensor<Real, 3> getU();
    
    protected:
        Eigen::Tensor<Real, 3> U_;
        Option<Real> option_;
};

template<std::floating_point Real>
class FixedStencil : public PDESolver<Real> {
    public:
        // Constructor
        FixedStencil();

        void _solve();

    private:

        double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2);
        void solve_sparse_matrix(std::vector<std::vector<double>>& U);
        void policy_iteration(std::vector<std::vector<double>>& U);
        void apply_boundary_conditions(std::vector<std::vector<double>>& U);
};


#endif //HH__PDESOLVERSHPP_HH