#include "Option.hpp"

template<std::floating_point Real>
Option<Real>::Option(std::function<Real(Real, Real)> payoff,
          const Real r, const std::pair<Real, Real>& sigmas_1,
          const std::pair<Real, Real>& sigmas_2, 
          const std::pair<Real, Real>& rho, const Real expiry)
:payoff_(payoff), r_(r), sigmas_1_(sigmas_1), sigmas_2_(sigmas_2),
rho_(rho), expiry_(expiry) {}

template<std::floating_point Real>
Real Option<Real>::evaluate(const Real S1, const Real S2)
{
    return payoff_(S1, S2);
}

template<std::floating_point Real>
Real Option<Real>::getExpiry()
{
    return expiry_;
}

template<std::floating_point Real>
Real Option<Real>::getDiscountRate()
{
    return r_;
}

template<std::floating_point Real>
std::pair<Real, Real> Option<Real>::getSigmas_1()
{
    return sigmas_1_;
}

template<std::floating_point Real>
std::pair<Real, Real> Option<Real>::getSigmas_2()
{
    return sigmas_2_;
}

template<std::floating_point Real>
std::pair<Real, Real> Option<Real>::getRhos()
{
    return rho_;
}

template class Option<double>;