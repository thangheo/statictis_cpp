#include <iostream>
#include <vector>
#include <unordered_map>
#include <cmath>

double calculateNormalizedEntropy(const std::vector<double>& data) {
    int dataSize = data.size();
    std::unordered_map<double, int> frequencyMap;
    
    // Count the frequency of each unique value in the time series
    for (const double& value : data) {
        frequencyMap[value]++;
    }
    
    double normalizedEntropy = 0.0;
    
    // Calculate the probability of each unique value and the corresponding entropy contribution
    for (const auto& pair : frequencyMap) {
        double probability = static_cast<double>(pair.second) / dataSize;
        normalizedEntropy -= probability * log2(probability);
    }
    
    // Normalize the entropy between 0 and 1
    normalizedEntropy /= log2(dataSize);
    return normalizedEntropy;
}

int main() {
    // Create a vector of time series data
    std::vector<double> data = { -0,1,2,-2,4,9,7,6 ,-0,1,2,-2,4,9,7,6,-0,1,2,-2,4,9,7,6};

    // Calculate the normalized entropy
    double normalizedEntropy = calculateNormalizedEntropy(data);
    
    // Output the result
    std::cout << "Normalized Entropy: " << normalizedEntropy << std::endl;
    
    return 0;
}
