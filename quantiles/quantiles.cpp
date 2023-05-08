#include <iostream>
#include <vector>
#include <algorithm>

int main() {
    // Create a vector of data
    std::vector<double> data = { 1.2, 2.4, 3.6, 4.8, 6.0 };

    // Calculate quartiles
    std::sort(data.begin(), data.end());
    double firstQuartile = data[data.size() * 0.25];
    double secondQuartile = data[data.size() * 0.50];
    double thirdQuartile = data[data.size() * 0.75];

    // Output the results
    std::cout << "First Quartile: " << firstQuartile << std::endl;
    std::cout << "Second Quartile (Median): " << secondQuartile << std::endl;
    std::cout << "Third Quartile: " << thirdQuartile << std::endl;

    return 0;
}
