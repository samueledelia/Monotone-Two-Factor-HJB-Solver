# Monotone Two-Factor HJB Solver

This project implements a fully implicit, unconditionally monotone numerical scheme to solve the two-factor uncertain volatility model. The solver is designed to ensure stability and convergence to the viscosity solution. It focuses on pricing financial derivatives under worst-case volatility scenarios, employing a hybrid discretization method to achieve robust and accurate results.

You can read the original paper from [here](https://cs.uwaterloo.ca/~paforsyt/uncertain_2d.pdf)

## Features

- **Unconditionally Monotone Scheme**: Ensures numerical stability and adherence to the theoretical properties of the model.
- **Hybrid Discretization**: Combines finite difference methods to handle the two-factor uncertainty efficiently.
- **Financial Derivative Pricing**: Targets scenarios like European call options and butterfly options under uncertain volatility.

## Requirements

To run this project, you will need the following:

1. **Eigen Library**: A high-performance C++ library for linear algebra. You can find it [here](https://eigen.tuxfamily.org/).
2. **GoogleTest**: For unit testing the codebase. Clone it using:
```bash
   git clone git@github.com:google/googletest.git
```

### Cloning the Repository

Clone the repository from [here](https://github.com/diantonioandrea/pacs-project):

```bash
git clone git@github.com:diantonioandrea/pacs-project.git
```

### Compilation and Execution

The code is written to work with the C++ standard library.

Compile everything[^compilation] with:

[^compilation]: The `-j` option is advised.

```bash
make
```

Compile examples with:

```bash
make examples
```

Compile tests with:

```bash
make tests
```

## Using the Code

To use the code, include `<PDESolver.hpp>` which provides all necessary tools.



### Examples

Examples are divided into the following categories based on different option to price:

1. European Call Option on the maximun of two underlying asset: `max_eu_call.cpp`

2. Butterfly Option: `butterfly.cpp`