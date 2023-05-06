#ifndef MS_SHAPIRO_WILK_HPP_
#define MS_SHAPIRO_WILK_HPP_

#include <deque>
#include <utility>

class ms_shapiro_wilk
{
public:
    ms_shapiro_wilk();
    ms_shapiro_wilk(std::deque<std::pair<size_t, double>> x, long n, long n1, long n2);
    ms_shapiro_wilk(std::deque<std::pair<size_t, double>> x, long n);

    ~ms_shapiro_wilk();
    double getResult() const;
    double getPValue() const;
    double getErrorCode() const;
    void appendSampleValue(double x);
    void clearSampleValues();
    void calculate(long n, long n1, long n2);
    double poly(const double *cc, int nord, double x);
    double ppnd7(double p, int &ifault) ;
    double alnorm(double x,bool upper) ;
    long sign(long x,long y);
    void swilk(bool init, std::deque<std::pair<size_t, double>> x, long n, long n1, long n2, std::deque<double>& a, double& w, double& pw, int& ifault);

private:
    std::deque<std::pair<size_t, double>> x_;
    double w_;
    double pw_;
    int ifault_;
    std::deque<double> a_;

};

#endif  // MS_SHAPIRO_WILK_HPP_
