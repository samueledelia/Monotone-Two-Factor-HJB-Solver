#include <iostream>
#include <vector>
#include <cmath>
#include <algorithm>
#include <eigen3/Eigen/Sparse>

// TEST CASE: EU Call Option on the maximum of two asset
// Market constants
const double r = 0.05; // Risk-free rate
const double q1 = 0.02; // Dividend yield for S1
const double q2 = 0.02; // Dividend yield for S2

// Range for the uncertain parameters
const double sigma1_min = 0.1, sigma1_max = 0.3;
const double sigma2_min = 0.1, sigma2_max = 0.3;
const double rho_min = -1.0, rho_max = 1.0;

// Define grid parameters
const int N1 = 91; // Number of grid points for S1
const int N2 = 91; // Number of grid points for S2
const double S1_max = 100.0;
const double S2_max = 100.0;
const double dS1 = S1_max / N1;
const double dS2 = S2_max / N2;
const double T = 1.0; // Time to maturity
const double dt = 0.01;
const int Nt = T / dt;

// Define the payoff function (e.g., European call option on the maximum)
double payoff(double S1, double S2) {
    double K = 50.0; // Strike price
    return std::max(S1 + S2 - K, 0.0);
}

// Define the differential operators
double LQ_U(double U, double sigma1, double sigma2, double rho, double S1, double S2) {
    double term1 = 0.5 * sigma1 * sigma1 * S1 * S1 * (U / (dS1 * dS1));
    double term2 = 0.5 * sigma2 * sigma2 * S2 * S2 * (U / (dS2 * dS2));
    double term3 = rho * sigma1 * sigma2 * S1 * S2 * (U / (dS1 * dS2));
    double term4 = (r - q1) * S1 * (U / dS1);
    double term5 = (r - q2) * S2 * (U / dS2);
    return term1 + term2 + term3 + term4 + term5 - r * U;
}

Eigen::SparseMatrix<double> vector_to_sparse(const std::vector<std::vector<double>>& data) {
  int rows = data.size();
  int cols = data[0].size();

  std::vector<Eigen::Triplet<double>> triplets;
  triplets.reserve(rows * cols); // Reserve space for potential non-zero elements

  for (int i = 0; i < rows; ++i) {
    for (int j = 0; j < cols; ++j) {
      if (data[i][j] != 0.0) { // Only add non-zero elements
        triplets.push_back(Eigen::Triplet<double>(i, j, data[i][j]));
      }
    }
  }

  Eigen::SparseMatrix<double> sparse_matrix(rows, cols);
  sparse_matrix.setFromTriplets(triplets.begin(), triplets.end());
  return sparse_matrix;
}

// This is a placeholder for a sparse matrix solver using Bi-CGSTAB method
// In practice, you would need a library such as Eigen or an implementation for sparse matrix operations.
void solve_sparse_matrix(std::vector<std::vector<double>>& U) {
    // Use Bi-CGSTAB or other iterative methods to solve the matrix equations
}

// Define policy iteration algorithm
void policy_iteration(std::vector<std::vector<double>>& U) {
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
void apply_boundary_conditions(std::vector<std::vector<double>>& U) {
    for (int i = 0; i < N1; ++i) {
        U[i][0] = 0; // Boundary at S2 = 0
        U[i][N2 - 1] = payoff(i * dS1, S2_max); // Boundary at S2 = S2_max
    }
    for (int j = 0; j < N2; ++j) {
        U[0][j] = 0; // Boundary at S1 = 0
        U[N1 - 1][j] = payoff(S1_max, j * dS2); // Boundary at S1 = S1_max
    }
}

int main() {
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

    return 0;
}
