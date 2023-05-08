#include <iostream>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

int main() {
    // Create a vector of data
    std::vector<double> data = { 1.2, 2.4, 3.6, 4.8, 6.0 };

    // Create an accumulator to calculate the skewness
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::skewness>> acc;

    // Add data points to the accumulator
    for (const auto& value : data) {
        acc(value);
    }

    // Calculate the skewness
    double skewness = boost::accumulators::skewness(acc);

    // Output the result
    std::cout << "Skewness: " << skewness << std::endl;

    return 0;
}
