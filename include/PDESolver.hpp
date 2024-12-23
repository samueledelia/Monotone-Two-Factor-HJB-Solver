#ifndef HH__PDESOLVERSHPP_HH
#define HH__PDESOLVERSHPP_HH
#include <vector>
#include <functional>

// ! Numerical HJB PDE Solver Class
class PDESolver {
    private:

    public:
        virtual void solve() const=0;
};

class FixedStencil : public PDESolver {
    private:

    public:
        // Constructor
        FixedStencil();

        double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2);
        void solve_sparse_matrix(std::vector<std::vector<double>>& U);
        void policy_iteration(std::vector<std::vector<double>>& U);
        void apply_boundary_conditions(std::vector<std::vector<double>>& U);
        void solve();
};


// Define the differential operators
double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2) {
    double term1 = 0.5 * sigma1 * sigma1 * S1 * S1 * (U / (dS1 * dS1));
    double term2 = 0.5 * sigma2 * sigma2 * S2 * S2 * (U / (dS2 * dS2));
    double term3 = rho * sigma1 * sigma2 * S1 * S2 * (U / (dS1 * dS2));
    double term4 = (r - q1) * S1 * (U / dS1);
    double term5 = (r - q2) * S2 * (U / dS2);
    return term1 + term2 + term3 + term4 + term5 - r * U;
}

void solve_sparse_matrix(std::vector<std::vector<double>>& U) {
    // Use Bi-CGSTAB or other iterative methods to solve the matrix equations
}

// Define policy iteration algorithm
void FixedStencil::policy_iteration(std::vector<std::vector<double>>& U) {
    for (int n = Nt - 1; n >= 0; --n) {
        for (int i = 1; i < N1 - 1; ++i) {
            for (int j = 1; j < N2 - 1; ++j) {
                // Discretize control set and maximize the value
                double max_val = -1e9;
                for (double sigma1 = sigma1_min; sigma1 <= sigma1_max; sigma1 += (sigma1_max - sigma1_min) / 10) {
                    for (double sigma2 = sigma2_min; sigma2 <= sigma2_max; sigma2 += (sigma2_max - sigma2_min) / 10) {
                        for (double rho = rho_min; rho <= rho_max; rho += (rho_max - rho_min) / 10) {
                            double val = LQ_U(U[i][j], sigma1, sigma2, rho, i * dS1, j * dS2);
                            if (val > max_val) {
                                max_val = val;
                            }
                        }
                    }
                }
                U[i][j] = U[i][j] + dt * max_val;
            }
        }
    }
}

// Apply boundary conditions
void FixedStencil::apply_boundary_conditions(std::vector<std::vector<double>>& U) {
    for (int i = 0; i < N1; ++i) {
        U[i][0] = 0; // Boundary at S2 = 0
        U[i][N2 - 1] = payoff(i * dS1, S2_max); // Boundary at S2 = S2_max
    }
    for (int j = 0; j < N2; ++j) {
        U[0][j] = 0; // Boundary at S1 = 0
        U[N1 - 1][j] = payoff(S1_max, j * dS2); // Boundary at S1 = S1_max
    }
}

void FixedStencil::solve() {
    // Initialize the grid and set boundary conditions
    std::vector<std::vector<double>> U(N1, std::vector<double>(N2, 0.0));
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            U[i][j] = payoff(i * dS1, j * dS2);
        }
    }
    apply_boundary_conditions(U);

    // Perform policy iteration
    policy_iteration(U);

    // Solve the sparse matrix
    solve_sparse_matrix(U);

    // Output the results (for example, save to a file or print to console)
    for (int i = 0; i < N1; ++i) {
        for (int j = 0; j < N2; ++j) {
            std::cout << U[i][j] << " ";
        }
        std::cout << std::endl;
    }
}


#endif //HH__PDESOLVERSHPP_HH