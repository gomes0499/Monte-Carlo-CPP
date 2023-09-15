# Monte Carlo Portfolio Optimization

## Introduction

This project uses Monte Carlo simulation to optimize a portfolio of assets. It calculates the daily returns of the assets, optimizes the weights of the assets in the portfolio to maximize the Sharpe ratio, and finally returns the optimal weights for the assets.

## Monte Carlo Simulation

Monte Carlo simulation is a statistical method used to model the probability of different outcomes in a process that cannot easily be predicted due to the intervention of random variables. It's a technique used to understand the impact of risk and uncertainty in prediction and forecasting models.

In the context of portfolio optimization, Monte Carlo simulation can be used to simulate different combinations of weights for the assets and calculate the corresponding Sharpe ratio.

## Sharpe Ratio

The Sharpe ratio is a measure to calculate the risk-adjusted performance of an investment. The Sharpe ratio is the average return earned in excess of the risk-free rate per unit of volatility or total risk.

Sharpe Ratio = (Rx – Rf) / StdDev Rx

where ( R_x ) is the portfolio return, ( R_f ) is the risk-free rate, and \( StdDev Rx ) is the volatility of the portfolio.

## How the Program Works

1. **Read CSV**: The program starts by reading a CSV file containing historical price data of the assets. The file should have one column for dates and one column for each asset.

2. **Calculate Daily Returns**: The program calculates the daily returns for each asset based on the historical prices.

3. **Portfolio Optimization**: Using Monte Carlo simulation, the program tests various combinations of weights for the assets in the portfolio and picks the combination that maximizes the Sharpe ratio.

4. **Results**: The program returns the optimal weights for the assets in the portfolio.

## Requirements

- C++
- Eigen Library

## Customizing the CSV File

To use your own CSV file for historical price data, follow these steps:

1. Place your CSV file in an accessible directory.

2. Open the `main.cpp` file and locate the following line:

    ```cpp
    std::map<std::string, Eigen::VectorXd> data = readCsv("../price.csv", matrixData);
    ```

3. Change the file path to the location where your new CSV file is saved. For example, if your file is named `my_prices.csv` and it's in the `Desktop` folder, the line should look something like this:

    ```cpp
    std::map<std::string, Eigen::VectorXd> data = readCsv("/path/to/your/folder/my_prices.csv", matrixData);
    ```

4. Save the changes and recompile the program.



## Compilation and Execution

```bash
g++ -o Monte_Carlo_CPP main.cpp
./Monte_Carlo_CPP
