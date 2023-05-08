#include <iostream>
#include <vector>

// Calculate cross-correlation between two time series at a given lag
double calculateCrossCorrelation(const std::vector<double>& series1, const std::vector<double>& series2, int lag) {
    int size1 = series1.size();
    int size2 = series2.size();

    // Check if the lag value is within the range of the time series
    if (lag >= size1 || lag >= size2) {
        std::cerr << "Error: Lag value exceeds time series size." << std::endl;
        return 0.0;
    }

    // Calculate the mean of the two time series
    double mean1 = 0.0;
    double mean2 = 0.0;

    for (int i = 0; i < size1; i++) {
        mean1 += series1[i];
    }
    mean1 /= size1;

    for (int i = 0; i < size2; i++) {
        mean2 += series2[i];
    }
    mean2 /= size2;

    // Calculate the cross-correlation
    double crossCorrelation = 0.0;

    for (int i = 0; i < size1 - lag; i++) {
        crossCorrelation += (series1[i] - mean1) * (series2[i + lag] - mean2);
    }

    crossCorrelation /= (size1 - lag);

    return crossCorrelation;
}

int main() {
    // Example usage
    std::vector<double> series1 = {1.0, 2.0, 3.0, 4.0, 5.0};
    std::vector<double> series2 = {2.0, 4.0, 6.0, 8.0, 10.0};

    int lag = 3;

    double correlation = calculateCrossCorrelation(series1, series2, lag);

    std::cout << "Cross-correlation at lag " << lag << ": " << correlation << std::endl;

    return 0;
}
