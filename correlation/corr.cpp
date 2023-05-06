#include <iostream>
#include <vector>
#include <cmath>

#include <boost/math/special_functions/beta.hpp>
#include <boost/math/distributions.hpp>
#include <boost/accumulators/accumulators.hpp>
#include <boost/accumulators/statistics/covariance.hpp>
#include <boost/accumulators/statistics/stats.hpp>
#include <boost/accumulators/statistics.hpp>
#include <boost/accumulators/statistics/variance.hpp>
using namespace boost::accumulators;
double kendall_tau(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();
    double concordant = 0.0;
    double discordant = 0.0;
    double ties_x = 0.0;
    double ties_y = 0.0;
    double tied_pairs = 0.0;

    for (size_t i = 0; i < n - 1; ++i) {
        for (size_t j = i + 1; j < n; ++j) {
            if (x[i] == x[j] && y[i] == y[j]) {
                ++tied_pairs;
            } else if (x[i] == x[j]) {
                ++ties_x;
                if (y[i] < y[j]) {
                    ++concordant;
                } else {
                    ++discordant;
                }
            } else if (y[i] == y[j]) {
                ++ties_y;
                if (x[i] < x[j]) {
                    ++concordant;
                } else {
                    ++discordant;
                }
            } else {
                if ((x[i] < x[j] && y[i] < y[j]) || (x[i] > x[j] && y[i] > y[j])) {
                    ++concordant;
                } else {
                    ++discordant;
                }
            }
        }
    }

    double denominator = std::sqrt((concordant + discordant + tied_pairs) *
                                   (concordant + discordant + ties_x) *
                                   (concordant + discordant + ties_y));
    double numerator = 2.0 * (concordant - discordant);
    double tau = numerator / denominator;
    return tau;
}

double pearson_corr(const std::vector<double>& x, const std::vector<double>& y) {
    size_t n = x.size();

    double sum_x = 0.0, sum_y = 0.0, sum_xy = 0.0;
    double sum_x_squared = 0.0, sum_y_squared = 0.0;

    for (size_t i = 0; i < n; ++i) {
        sum_x += x[i];
        sum_y += y[i];
        sum_xy += x[i] * y[i];
        sum_x_squared += x[i] * x[i];
        sum_y_squared += y[i] * y[i];
    }

    double numerator = (n * sum_xy) - (sum_x * sum_y);
    double denominator = std::sqrt((n * sum_x_squared - sum_x * sum_x) * (n * sum_y_squared - sum_y * sum_y));

    if (denominator == 0) {
        return 0.0;
    } else {
        return numerator / denominator;
    }
}


int main() {
    std::vector<double> x = {1, 2, 3, 4, -5};
    std::vector<double> y = {1, 2, 3, 4, 6};

    double tau = kendall_tau(x, y);
    double p = pearson_corr(x, y);

    std::cout << "Kendall's tau correlation coefficient: " << tau << std::endl;
    std::cout << "pearson_correlation coefficient: " << p << std::endl;

    return 0;
}
