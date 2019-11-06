#include "Statistics.h"
#include <iomanip>

namespace ZIRAN {
Quantiler::Quantiler(size_t num_bins)
    : num_bins(num_bins)
    , num_samples(0)
{
}

void Quantiler::insert(double sample)
{
    num_samples++;
    auto upper = std::upper_bound(marker_height.begin(), marker_height.end(), sample);

    if (marker_height.size() <= num_bins) {
        marker_height.insert(upper, sample);
        marker_position.emplace_back(marker_height.size());
        return;
    }
    assert(marker_height.size() == marker_position.size());
    assert(marker_height.size() == num_bins + 1);

    if (upper == marker_height.begin()) {
        *upper = sample;
    }
    else if (upper == marker_height.end()) {
        upper--;
        *upper = sample;
    }

    for (size_t k = upper - marker_height.begin(); k < marker_position.size(); k++) {
        marker_position[k]++;
    }

    for (size_t i = 1; i < num_bins; i++) {
        double desired_position = 1 + i * double(num_samples + 1) / num_bins;

        double qim1 = marker_height[i - 1];
        double qi = marker_height[i];
        double qip1 = marker_height[i + 1];

        size_t nim1 = marker_position[i - 1];
        size_t ni = marker_position[i];
        size_t nip1 = marker_position[i + 1];

        size_t diffm1 = ni - nim1;
        size_t diffp1 = nip1 - ni;

        if (desired_position > ni && diffp1 > 1) {
            double q_quadratic = qi + ((diffm1 + 1) * (qip1 - qi) / diffp1 + (diffp1 - 1) * (qi - qim1) / diffm1) / (nip1 - nim1);
            if (q_quadratic < qim1 || q_quadratic > qip1)
                marker_height[i] = qi + (qip1 - qi) / diffp1;
            else
                marker_height[i] = q_quadratic;
            marker_position[i]++;
        }
        else if (desired_position < ni && diffm1 > 1) {
            double q_quadratic = qi - ((diffm1 - 1) * (qip1 - qi) / diffp1 + (diffp1 + 1) * (qi - qim1) / (diffm1)) / (nip1 - nim1);
            if (q_quadratic < qim1 || q_quadratic > qip1)
                marker_height[i] = qi + (qim1 - qi) / (diffm1);
            else
                marker_height[i] = q_quadratic;
            marker_position[i]--;
        }
    }
}

size_t Quantiler::closetsBinToQuantile(double q)
{
    return std::min(size_t(q * num_bins), num_bins);
}

double Quantiler::quantile(double q)
{
    return marker_height[closetsBinToQuantile(q)];
}

void Quantiler::print(std::ostream& out)
{
    if (marker_height.size() <= num_bins) {
        out << "Not enough data\n";
        return;
    }
    for (size_t i = 0; i < marker_position.size(); i++) {
        double quantile = (100 * marker_position[i]) / double(num_samples);
        out << std::setprecision(1) << std::fixed << std::setw(4) << std::left << quantile << "%: ";
        out << std::setprecision(3) << std::right << std::scientific << marker_height[i] << "\n";
    }
}

void Quantiler::constructCumulativeHistogram(StdVector<double>& histogram, double start_quantile, double end_quantile)
{
    if (marker_height.size() <= num_bins) {
        return;
    }

    size_t width = histogram.size();
    double data_start = marker_height[closetsBinToQuantile(start_quantile)];
    double data_end = marker_height[closetsBinToQuantile(end_quantile)];
    double data_width = data_end - data_start;
    std::cout << '[' << data_start << ',' << data_end << ']' << '\n';
    for (size_t i = 0; i < width; i++) {
        double sample = (data_width * i) / width + data_start;
        auto upper = std::upper_bound(marker_height.begin(), marker_height.end(), sample);
        size_t k = upper - marker_height.begin();
        double alpha = (*upper - sample) / (*upper - *(upper - 1));
        histogram[i] = (alpha * marker_position[k - 1] + (1.0 - alpha) * marker_position[k]) / num_samples;
    }
}

/**
  Constructs an approximate histogram by using central differences on the cumulative one
  */
double Quantiler::constructHistogram(StdVector<double>& histogram, double start_quantile, double end_quantile)
{
    size_t width = histogram.size();
    constructCumulativeHistogram(histogram, start_quantile, end_quantile);
    double max = 0.0;
    for (size_t i = 0; i < width; i++) {
        double prev = start_quantile;
        double next = end_quantile;
        if (i > 0)
            prev = histogram[i - 1];
        if (i + 1 < width)
            next = histogram[i + 1];
        double density = 0.5 * (next - prev) * width;
        if (density > max)
            max = density;
        histogram[i] = density;
    }
    return max;
}

void Quantiler::printGraph(std::ostream& out, const StdVector<double>& y, double data_height, size_t plot_width, size_t plot_height)
{
    double scale = plot_height / data_height;
    for (size_t i = 0; i < plot_height; i++) {
        for (size_t j = 0; j < plot_width; j++)
            out << ((scale * y[j] < (plot_height - i)) ? ' ' : '*');
        out << '\n';
    }
}

void Quantiler::printScale(std::ostream& out, size_t width, double start_quantile, double end_quantile)
{
    size_t field_width = 16;
    out << std::left << std::setprecision(1) << std::scientific;
    size_t d = width / field_width;
    for (size_t i = 0; i < d - 1; i++) {
        size_t k = closetsBinToQuantile(start_quantile + (end_quantile - start_quantile) * i / d);
        out << std::setw(field_width) << marker_height[k];
    }
    out << std::right << std::setw(field_width) << marker_height[closetsBinToQuantile(end_quantile)];
    out << std::endl;
}

void Quantiler::printCumulativeHistogram(std::ostream& out, size_t plot_width, size_t plot_height, double start_quantile, double end_quantile)
{
    if (marker_height.size() <= num_bins) {
        out << "Not enough data\n";
        return;
    }

    StdVector<double> histogram(plot_width);
    constructCumulativeHistogram(histogram, start_quantile, end_quantile);
    printGraph(out, histogram, end_quantile, plot_width, plot_height);
    printScale(out, plot_width, start_quantile, end_quantile);
}

void Quantiler::printHistogram(std::ostream& out, size_t plot_width, size_t plot_height, double start_quantile, double end_quantile)
{
    if (marker_height.size() <= num_bins) {
        out << "Not enough data\n";
        return;
    }

    StdVector<double> histogram(plot_width);
    double max = constructHistogram(histogram, start_quantile, end_quantile);
    printGraph(out, histogram, max, plot_width, plot_height);
    printScale(out, plot_width, start_quantile, end_quantile);
}
} // namespace ZIRAN
