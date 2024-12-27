#include "PDESolver.hpp"
#include <eigen3/Eigen/Dense>
#include <eigen3/Eigen/Sparse>
#include <eigen3/Eigen/SparseLU>

constexpr uint16_t BS_T_DIM = 2;

using Vector = Eigen::VectorXd;

template<std::floating_point Real>
class BSSolver : public PDESolver<Real, BS_T_DIM>
{
public:
    BSSolver(uint32_t N1, uint32_t N_tau,
             std::unique_ptr<OneAssetOption<Real>> option, Real S_max);

    void solve();

    Eigen::Tensor<Real, BS_T_DIM>& getU() override;

    Eigen::MatrixXd getV();

protected:
    Eigen::Tensor<Real, BS_T_DIM> initBoundaryConditions();
    Eigen::MatrixXd initBoundaryConditions_();
    Real getSigma();

    const uint32_t N_;
    const double dS_, dt_;
    const Real S_max_;
    Vector S_;
    Eigen::MatrixXd V_;
private:
    Eigen::Tensor<Real, BS_T_DIM> MatrixToTensor(const Eigen::MatrixXd& matrix);
};