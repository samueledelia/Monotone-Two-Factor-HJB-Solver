
#ifndef OPTIONBASE_H
#define OPTIONBASE_H
#include <concepts>
#include <functional>

template<std::floating_point Real>
class OptionBase {
public:
    OptionBase(Real r, Real expiry)
        :r_(r), expiry_(expiry) {}

    Real getExpiry() { return expiry_; }

    Real getDiscountRate() { return r_; }

    template<typename... Reals>
    Real evaluate(Reals... reals)
    {
        return evaluate_(std::index_sequence_for<Reals...>(), std::forward<Reals>(reals)...);
    }

    //FIXME: temporary, better solution has to be found
    virtual Real getSigma() { return 0.0; }

    virtual std::pair<Real, Real> getSigmas_1() { return std::make_pair(0.0, 0.0); }

    virtual std::pair<Real, Real> getSigmas_2() { return std::make_pair(0.0, 0.0); }

private:

    template <std::size_t... Is, typename... Args>
    Real evaluate_(std::index_sequence<Is...>, Args&&... args) {
        // Compile-time branching based on the number of arguments
        if constexpr (sizeof...(Args) == 0) {
            static_assert(sizeof...(Args) == 0, "No arguments provided.");
            return -1; // No arguments
        } else if constexpr (sizeof...(Args) == 1) {
            return evaluateOption_(std::forward<Args>(args)...); // One argument
        } else if constexpr (sizeof...(Args) == 2) {
            return evaluateOption_(std::forward<Args>(args)...); // Two arguments
        } else {
            // Handle the case where the number of arguments is not supported
            static_assert(sizeof...(Args) <= 3, "Too many arguments provided."); // Compile time error
            return -1; // Or throw an exception at runtime if needed.
        }
    }

protected:
    virtual Real evaluateOption_([[maybe_unused]] Real a) {
        return 0.0;
    }
    virtual Real evaluateOption_([[maybe_unused]] Real a,[[maybe_unused]]  Real b) {
        return 0.0;
    }

    Real r_;
    Real expiry_;
};

template class OptionBase<double>;

#endif //OPTIONBASE_H
