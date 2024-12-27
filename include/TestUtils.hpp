#include <fstream>
#include <sstream>
#include <vector>
#include <eigen3/Eigen/Dense>
#include <gtest/gtest.h>

template<typename MType>
MType readMatrixFromFile(const std::string& filename) {
    std::ifstream file(filename);
    if (!file.is_open()) {
        throw std::runtime_error("Could not open the file!");
    }

    std::vector<std::vector<double>> matrix_data;
    std::string line;

    while (std::getline(file, line)) {
        std::stringstream ss(line);
        std::vector<double> row;
        double value;

        while (ss >> value) {
            row.push_back(value);
        }

        if (!row.empty()) {
            matrix_data.push_back(row);
        }
    }

    file.close();

    size_t rows = matrix_data.size();
    size_t cols = matrix_data.empty() ? 0 : matrix_data[0].size();

    MType matrix(rows, cols);

    for (size_t i = 0; i < rows; ++i) {
        for (size_t j = 0; j < cols; ++j) {
            matrix(i, j) = matrix_data[i][j];
        }
    }

    return matrix;
}

template<template <typename, int> class TensorType, typename Real, int Dim>
testing::AssertionResult TensorAreApproxEqual(
    const TensorType<Real, Dim>& expected,
    const TensorType<Real, Dim>& actual,
    const Real tolerance = 1e-5)
{
     for (int i = 0; i < expected.size(); ++i) {
        if (std::abs(expected.data()[i] - actual.data()[i]) > tolerance) {
            return testing::AssertionFailure()
                   << "Mismatch at index " << i
                   << ": expected " << expected.data()[i]
                   << " but got " << actual.data()[i];
        }
    }
    return testing::AssertionSuccess();
}
