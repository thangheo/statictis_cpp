#include <iostream>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>

using namespace boost::accumulators;

double calculateCovariance(const std::vector<double>& x, const std::vector<double>& y) {
    accumulator_set<double, stats<tag::covariance<double, tag::covariate1> > > acc;
    
    // Check if the sizes of the two datasets match
    if (x.size() != y.size()) {
        throw std::runtime_error("Datasets have different sizes");
    }
    
    // Iterate over the datasets and update the accumulator
    for (std::size_t i = 0; i < x.size(); ++i) {
        acc(x[i], covariate1 = y[i]);
    }
    
    // Calculate and return the covariance
    return covariance(acc);
}

int main() {
    // Example datasets
    std::vector<double> x = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> y = {2.0, 4.0, 6.0, 8.0, 10.0};
    
    try {
        double covariance = calculateCovariance(x, y);
        std::cout << "Covariance: " << covariance << std::endl;
    } catch (const std::runtime_error& e) {
        std::cout << "Error: " << e.what() << std::endl;
    }
    
    return 0;
}

