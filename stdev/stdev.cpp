#include <iostream>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/math/statistics/univariate_statistics.hpp>
// int main() {
//     // Create a vector of data
//     std::vector<double> data = { 1.2, 2.4, 3.6, 4.8, 6.0 };

//     // Calculate the standard deviation using Boost
//     double standardDeviation = boost::math::statistics::standard_deviation(data.begin(), data.end());

//     // Output the result
//     std::cout << "Standard Deviation: " << standardDeviation << std::endl;

//     return 0;
// }


int main() {
    // Create a vector of data
    std::vector<double> data = { 1.2, 2.4, 3.6, 4.8, 6.0 };

    // Create an accumulator to calculate the standard deviation
    boost::accumulators::accumulator_set<double, boost::accumulators::stats<boost::accumulators::tag::variance>> acc;

    // Add data points to the accumulator
    for (const auto& value : data) {
        acc(value);
    }

    // Calculate the standard deviation
    double standardDeviation = std::sqrt(boost::accumulators::variance(acc));

    // Output the result
    std::cout << "Standard Deviation: " << standardDeviation << std::endl;

    return 0;
}
