#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <map>
#include <vector>
#include <random>
#include <cmath>
#include <Eigen/Dense>

Eigen::MatrixXd mapToEigenMatrix(const std::map<std::string, Eigen::VectorXd> &data) {
    // Assuming all vectors have the same length
    int rows = data.begin()->second.size();
    int cols = data.size();

    Eigen::MatrixXd matrix(rows, cols);

    int col = 0;
    for (const auto &pair: data) {
        matrix.col(col) = pair.second;
        col++;
    }

    return matrix;
}

// Function to read the CSV.
std::map<std::string, Eigen::VectorXd> readCsv(const std::string &filename, Eigen::MatrixXd &matrixData) {
    std::map<std::string, Eigen::VectorXd> data;
    std::ifstream file(filename);
    std::string line;

    if (!file.is_open()) {
        std::cout << "Error: Could not open file " << filename << std::endl;
        return data;
    }

    // Read the header to get the stock symbols
    std::getline(file, line);
    std::stringstream headerStream(line);
    std::string cell;
    std::vector<std::string> symbols;

    // Skip date column
    std::getline(headerStream, cell, ',');

    while (std::getline(headerStream, cell, ',')) {
        symbols.push_back(cell);
        data[cell] = Eigen::VectorXd();  // Initialize a VectorXd for each symbol
    }

    std::vector<std::vector<double>> tempData(symbols.size());

    // Read the rest of the lines
    while (std::getline(file, line)) {
        std::stringstream lineStream(line);

        // Skip date column for each row
        std::getline(lineStream, cell, ',');

        for (size_t i = 0; i < symbols.size(); ++i) {
            std::getline(lineStream, cell, ',');
            tempData[i].push_back(std::stod(cell));
        }
    }

// Convert vector of vectors to Eigen::VectorXd and store it in the map
    for (size_t i = 0; i < symbols.size(); ++i) {
        Eigen::VectorXd vec(tempData[i].size());
        for (size_t j = 0; j < tempData[i].size(); ++j) {
            vec(j) = tempData[i][j];
        }
        data[symbols[i]] = vec;
    }

    matrixData = mapToEigenMatrix(data);  // populate matrixData
    return data;
}

// Function to calculate daily return
Eigen::VectorXd calculateDailyReturn(const Eigen::VectorXd &prices) {
    int n = prices.size() - 1;  // Number of daily returns
    Eigen::VectorXd dailyReturn(n);

    for (int i = 0; i < n; ++i) {
        dailyReturn(i) = (prices(i + 1) - prices(i)) / prices(i);
    }

    return dailyReturn;
}

// Function to calculate Standard Deviation
double calculateStandardDeviation(const Eigen::VectorXd &data) {
    int n = data.size();

    // If there's not enough data to calculate standard deviation, return 0
    if (n <= 1) {
        return 0.0;
    }

    // Calculate the mean of the data
    double mean = data.mean();

    // Calculate the sum of squared differences from the mean
    double sumOfSquares = 0.0;
    for (int i = 0; i < n; ++i) {
        sumOfSquares += std::pow(data(i) - mean, 2);
    }

    // Calculate the standard deviation
    double standardDeviation = std::sqrt(sumOfSquares / (n - 1));

    return standardDeviation;
}

// Function to calculate Monte Carlo Simulation
Eigen::VectorXd monteCarloSimulation(const Eigen::VectorXd &historicalPrices, int numSimulations) {
    int numDays = historicalPrices.size();

    // Calculate daily returns
    Eigen::VectorXd dailyReturns(numDays - 1);
    for (int i = 0; i < numDays - 1; ++i) {
        dailyReturns(i) = (historicalPrices(i + 1) - historicalPrices(i)) / historicalPrices(i);
    }

    // Calculate the mean and standard deviation of daily returns
    double meanReturn = dailyReturns.mean();
    double stdDevReturn = std::sqrt((dailyReturns.array() - meanReturn).square().sum() / (numDays - 2));

    // Initialize random number generation
    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    // Initialize the simulated prices vector
    Eigen::VectorXd simulatedPrices(numSimulations);

    // Monte Carlo simulation
    for (int i = 0; i < numSimulations; ++i) {
        double lastPrice = historicalPrices(numDays - 1);
        double z = distribution(generator);
        double nextPrice = lastPrice * std::exp(
                (meanReturn - 0.5 * stdDevReturn * stdDevReturn) + stdDevReturn * std::sqrt(1.0) * z);
        simulatedPrices(i) = nextPrice;
    }

    return simulatedPrices;
}

// Function to optimize the portfolio
Eigen::VectorXd optimizePortfolio(const Eigen::MatrixXd &returns) {
    const int NUM_SIMULATIONS = 1000000;
    const double RISK_FREE_RATE = 0.01;  // Assume a risk-free rate of 1%

    Eigen::VectorXd best_weights(3);
    double max_sharpe_ratio = -std::numeric_limits<double>::infinity();

    std::default_random_engine generator;

    for (int i = 0; i < NUM_SIMULATIONS; ++i) {
        Eigen::VectorXd weights(3);

        // Generate random weights that sum to 1
        double w1 = ((double) rand() / (RAND_MAX));
        double w2 = ((double) rand() / (RAND_MAX)) * (1 - w1);
        double w3 = 1 - w1 - w2;
        weights << w1, w2, w3;

        Eigen::MatrixXd portfolioReturn = returns * weights;

        double expected_return = portfolioReturn.mean();
        double risk = std::sqrt((portfolioReturn.array() - expected_return).square().mean());

        double sharpe_ratio = (expected_return - RISK_FREE_RATE) / risk;

        if (sharpe_ratio > max_sharpe_ratio) {
            max_sharpe_ratio = sharpe_ratio;
            best_weights = weights;
        }
    }

    return best_weights;
}


int main() {
    // Step 1: Read the CSV stocks data
    Eigen::MatrixXd matrixData;
    std::map<std::string, Eigen::VectorXd> data = readCsv("/Users/gomes/Desktop/Projects/Monte-Carlo-CPP/price.csv",
                                                          matrixData);
    std::cout << "Rows: " << matrixData.rows() << ", Columns: " << matrixData.cols() << std::endl;

    // Step 2: Calculate daily returns
    Eigen::MatrixXd returns = calculateDailyReturn(matrixData);

    // Step 3: Optimize portfolio
    Eigen::VectorXd best_weights = optimizePortfolio(returns);

    // Display the results
    std::cout << "Optimized Portfolio Weights:" << std::endl;
    std::cout << "Asset 1: " << best_weights(0) << std::endl;
    std::cout << "Asset 2: " << best_weights(1) << std::endl;
    std::cout << "Asset 3: " << best_weights(2) << std::endl;

    return 0;
}

