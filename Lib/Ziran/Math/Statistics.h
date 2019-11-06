#ifndef STATISTICS_H
#define STATISTICS_H

#include <vector>
#include <iostream>
#include <Ziran/CS/Util/Forward.h>

namespace ZIRAN {

/**
  Computes estimates for equal probability quantiles of a stream of data
  Used the P squared algorithm of Raj Jain and Iimrich Chlamtac

*/
class Quantiler {
public:
    size_t num_bins;
    size_t num_samples;
    StdVector<double> marker_height; // q
    StdVector<size_t> marker_position; // n

    Quantiler(size_t num_bins);

    void insert(double sample);

    size_t closetsBinToQuantile(double q);

    double quantile(double q);

    void print(std::ostream& out);

    void constructCumulativeHistogram(StdVector<double>& histogram, double start_quantile = 0.0, double end_quantile = 1.0);

    /**
      Constructs an approximate histogram by using central differences on the cumulative one
      */
    double constructHistogram(StdVector<double>& histogram, double start_quantile = 0.0, double end_quantile = 1.0);

    void printGraph(std::ostream& out, const StdVector<double>& y, double data_height, size_t plot_width, size_t plot_height);

    void printScale(std::ostream& out, size_t width, double start_quantile = 0.0, double end_quantile = 1.0);

    void printCumulativeHistogram(std::ostream& out, size_t plot_width, size_t plot_height, double start_quantile = 0.0, double end_quantile = 1.0);

    void printHistogram(std::ostream& out, size_t plot_width, size_t plot_height, double start_quantile = 0.0, double end_quantile = 1.0);
};
} // namespace ZIRAN
#endif
