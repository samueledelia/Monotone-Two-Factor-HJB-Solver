#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <vector>
#include <functional>
#include <Eigen/Sparse>


// ! Numerical HJB PDE Solver Class
class PDESolver {
    private:

    public:
        virtual void solve() const=0;
};

class FixedStencil : public PDESolver {
    public:
        // Constructor
        FixedStencil();

        void solve();

    private:

        double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2);
        void solve_sparse_matrix(std::vector<std::vector<double>>& U);
        void policy_iteration(std::vector<std::vector<double>>& U);
        void apply_boundary_conditions(std::vector<std::vector<double>>& U);
};


#endif //HH__PDESOLVERSHPP_HH