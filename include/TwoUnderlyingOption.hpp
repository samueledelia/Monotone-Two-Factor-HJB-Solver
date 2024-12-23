#ifndef HH__TWOUNDERLYINGOPTION_HH
#define HH__TWOUNDERLYINGOPTION_HH

class TwoUnderlyingOption {
public:
    virtual double payoff(double price1, double price2) const = 0; 
};



#endif // HH__TWOUNDERLYINGOPTION_HH