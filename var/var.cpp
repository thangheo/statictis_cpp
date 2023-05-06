#include <iostream>
#include <vector>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics.hpp>
#include "var.h"
int main() {
    std::vector<float> data {1.0, 2.0, 3.0, 4.0, 5.0,1000};
    std::cout << "Variance: " << variance(data) << std::endl;

    return 0;
}
float variance(std::vector<float> &data){
    boost::accumulators::accumulator_set<float, boost::accumulators::stats<boost::accumulators::tag::variance>> acc;

    for (auto& x : data) {
        acc(x);
    }
    return boost::accumulators::variance(acc) ;
}