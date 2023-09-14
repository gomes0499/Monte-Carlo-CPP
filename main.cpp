#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <map>
#include <cmath>
#include <random>
#include <Eigen/Dense>


// Function to read CSV
std::map<std::string, Eigen::VectorXd> readCsv(const std::string &filename) {
    std::ifstream file(filename);
    std::string line, word;
    std::vector<std::string> headers;
    std::map<std::string, std::vector<double>> tempData; // temporary storage for data
    std::map<std::string, Eigen::VectorXd> data;

    if (file.good()) {
        std::getline(file, line);
        std::stringstream ss(line);

        // Read Headers
        while (std::getline(ss, word, ',')) {
            headers.push_back(word);
        }

        // Read data
        while (std::getline(file, line)) {
            std::stringstream ss(line);
            std::getline(ss, word, ','); // Ignoring Date

            for (size_t i = 1; i < headers.size(); ++i) {
                std::getline(ss, word, ',');
                tempData[headers[i]].push_back(std::stod(word));
            }
        }

        // Convert std::vector<double> to Eigen::VectorXd
        for (const auto &item : tempData) {
            Eigen::VectorXd vec(item.second.size());
            vec = Eigen::VectorXd::Map(item.second.data(), item.second.size());
            data[item.first] = vec;
        }
    } else {
        std::cout << "Error opening the file." << std::endl;
    }

    return data;
}

// Function to calculate de daily return
Eigen::VectorXd calculateDailyReturn(Eigen::Matrix<double, -1, 1> prices) {
    Eigen::VectorXd returns(prices.size() - 1);
    for (size_t i = 1; i < prices.size(); ++i) {
        double dailyReturn = (prices[i] - prices[i - 1]) / prices[i - 1];
        returns(i - 1) = dailyReturn;
    }
    return returns;
}

// Function to calculate mean
double mean(const Eigen::VectorXd &v) {
    return v.mean();
}

// Function to calculate Standard Deviation
double std_dev(const Eigen::VectorXd &v, double mean) {
    return std::sqrt((v.array() - mean).square().mean());
}

// Function to calculate the Monte Carlo Simulation
Eigen::VectorXd monteCarloSimulation(const Eigen::VectorXd &prices, int numSimulations) {
    Eigen::VectorXd dailyReturn = calculateDailyReturn(prices);
    double meanReturn = dailyReturn.mean();
    double volatility = std::sqrt((dailyReturn.array() - meanReturn).square().mean());

    Eigen::VectorXd simulatedPrices(numSimulations);
    double currentPrice = prices(prices.size() - 1);

    std::default_random_engine generator;
    std::normal_distribution<double> distribution(0.0, 1.0);

    for (int i = 0; i < numSimulations; ++i) {
        double z = distribution(generator);
        double deltaT = 1;
        double drift = meanReturn - (std::pow(volatility, 2) / 2) * deltaT;
        double change = std::exp(drift * deltaT + volatility * std::sqrt(deltaT) * z);
        double newPrice = currentPrice * change;
        simulatedPrices(i) = newPrice;
    }
    return simulatedPrices;
}

// Function to optimize portfolio
Eigen::MatrixXd simulatePortfolio(const Eigen::MatrixXd &returns, const Eigen::VectorXd &weights) {
    std::cout << "returns: " << returns.rows() << " x " << returns.cols() << std::endl;
    std::cout << "weights: " << weights.rows() << " x " << weights.cols() << std::endl;

    Eigen::MatrixXd portfolioReturn = returns.transpose() * weights;
    Eigen::MatrixXd portfolioRisk = weights.transpose() * returns * weights;

    Eigen::MatrixXd result(2, 1);
    result(0, 0) = portfolioReturn.mean();
    result(1, 0) = std::sqrt(portfolioRisk.sum());

    return result;
}

int main() {
    // Lê os dados dos preços a partir do CSV
    auto data = readCsv("/Users/gomes/Desktop/Projects/Monte-Carlo-CPP/price.csv");
    std::map<std::string, Eigen::VectorXd> simulated_prices;  // Assumindo que monteCarloSimuation agora retorna Eigen::VectorXd
    const int NUM_SIMULATIONS = 1000;

    // Calcula as simulações de Monte Carlo para cada ativo e armazena em 'simulated_prices'
    for (const auto &asset: data) {
        simulated_prices[asset.first] = monteCarloSimulation(asset.second, NUM_SIMULATIONS);
    }


    // Cria uma matriz para armazenar os retornos simulados
    Eigen::MatrixXd returns(simulated_prices.size(), NUM_SIMULATIONS);

    // Popula a matriz 'returns' com os retornos simulados
    int i = 0;
    for (const auto &asset: simulated_prices) {
        Eigen::VectorXd dailyReturn = calculateDailyReturn(asset.second);
        returns.row(i) = dailyReturn;
        ++i;
    }

    double max_return = -std::numeric_limits<double>::infinity();
    Eigen::VectorXd best_weights(3);

    // Gerar carteiras aleatórias
    for (int i = 0; i < NUM_SIMULATIONS; ++i) {
        Eigen::VectorXd weights(3);

        // Gerar pesos aleatórios que somam 1
        double w1 = ((double) rand() / (RAND_MAX));
        double w2 = ((double) rand() / (RAND_MAX)) * (1 - w1);
        double w3 = 1 - w1 - w2;
        weights << w1, w2, w3;

        Eigen::MatrixXd result = simulatePortfolio(returns, weights);
        if (result(0, 0) > max_return) {
            max_return = result(0, 0);
            best_weights = weights;
        }
    }

    std::cout << "Melhor Carteira:" << std::endl;
    std::cout << "Peso do Ativo 1: " << best_weights(0) << std::endl;
    std::cout << "Peso do Ativo 2: " << best_weights(1) << std::endl;
    std::cout << "Peso do Ativo 3: " << best_weights(2) << std::endl;
    std::cout << "Retorno Esperado: " << max_return << std::endl;

    return 0;
}
