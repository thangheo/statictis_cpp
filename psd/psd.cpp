#include <iostream>
#include <vector>
#include <complex>
#include <cmath>
#include <algorithm>  // Include the <algorithm> header

// Apply a windowing function to the time series data
void applyWindow(std::vector<double>& data) {
    int dataSize = data.size();
    for (int i = 0; i < dataSize; i++) {
        double windowValue = 0.54 - 0.46 * cos(2 * M_PI * i / (dataSize - 1)); // Hamming window
        data[i] *= windowValue;
    }
}

// Calculate the Power Spectral Density (PSD) using FFT
std::vector<double> calculatePSD(const std::vector<double>& data) {
    int dataSize = data.size();

    // Apply windowing function to the time series data
    std::vector<double> windowedData = data;
    applyWindow(windowedData);

    // Perform the Fast Fourier Transform (FFT)
    std::vector<std::complex<double>> fftData(dataSize);
    for (int i = 0; i < dataSize; i++) {
        fftData[i] = std::complex<double>(windowedData[i], 0.0);
    }
    std::vector<std::complex<double>> fftResult(dataSize);
    std::transform(fftData.begin(), fftData.end(), fftResult.begin(), [&](const std::complex<double>& value) {
        // return std::abs(std::fft(value));
            return std::abs(value);
    });

    // Calculate the power spectrum (squared magnitude)
    std::vector<double> powerSpectrum(dataSize);
    std::transform(fftResult.begin(), fftResult.end(), powerSpectrum.begin(), [&](const std::complex<double>& value) {
        double magnitude = std::abs(value);
        return magnitude * magnitude;
    });

    // Normalize the power spectrum
    double normalizationFactor = 1.0 / (dataSize * dataSize);
    std::transform(powerSpectrum.begin(), powerSpectrum.end(), powerSpectrum.begin(), [&](double value) {
        return value * normalizationFactor;
    });

    return powerSpectrum;
}

int main() {
    // Create a vector of time series data
    std::vector<double> data = { 1.1,2,3,4,5,6,7,8,9,8,7,6,2,1.1,2,3,4,5,6,7,8,9,8,7,6,2,1.1,2,3,4,5,6,7,8,9,8,7,6,2 };

    // Calculate the Power Spectral Density (PSD)
    std::vector<double> psd = calculatePSD(data);

    // Output the PSD values
    for (int i = 0; i < psd.size(); i++) {
        double frequency = static_cast<double>(i) / data.size();
        std::cout << "Frequency: " << frequency << ", Power Spectral Density: " << psd[i] << std::endl;
    }

    return 0;
}
