#include <iostream>
#include <vector>
#include <cmath>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
using namespace std;
double calculateHurstExponent(const std::vector<double>& data) {
    int dataSize = data.size();
    double hurstExponent = 0.0;

    for (int blockSize = 2; blockSize <= dataSize; blockSize++) {
        int numBlocks = dataSize / blockSize;
        double avgRange = 0.0;
        double avgStdDev = 0.0;

        for (int i = 0; i < numBlocks; i++) {
            int startIndex = i * blockSize;
            int endIndex = startIndex + blockSize;
            std::vector<double> block(data.begin() + startIndex, data.begin() + endIndex);

            double blockRange = *std::max_element(block.begin(), block.end()) - *std::min_element(block.begin(), block.end());
            double blockStdDev = std::sqrt(std::accumulate(block.begin(), block.end(), 0.0, [](double sum, double value) {
                return sum + (value * value);
            }) / blockSize);

            avgRange += blockRange;
            avgStdDev += blockStdDev;
            
        }

        avgRange /= numBlocks;
        avgStdDev /= numBlocks;
        cout << "avgRange: " << avgRange << endl;
        cout << "avgStdDev: " << avgStdDev << endl;

        hurstExponent += std::log(avgRange / avgStdDev) / std::log(blockSize);
    }

    hurstExponent /= dataSize;

    return hurstExponent;
}


int main() {
    // Create a vector of data
    std::vector<double> data = { 1.1,2,3,4,5,6,7,8,9,8,7,6,2,1.1,2,3,4,5,6,7,8,9,8,7,6,2,1.1,2,3,4,5,6,7,8,9,8,7,6,2 };

    // Calculate the Hurst Exponent
    double hurstExponent = calculateHurstExponent(data);

    // Output the result
    std::cout << "Hurst Exponent: " << hurstExponent << std::endl;

    return 0;
}
