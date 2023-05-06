#include <cmath>
#include <limits>
#include <algorithm>
#include <iostream>
#include "shapiro.h"
#include <numeric> 
//the original source I take from this page : http://www.matrixscience.com/msparser/help/group___shapiro_wilk_source_code.html


ms_shapiro_wilk::ms_shapiro_wilk()
    :w_(0)
    ,pw_(0)
    ,ifault_(0) 
{
}

ms_shapiro_wilk::ms_shapiro_wilk(std::deque<std::pair<size_t,double> > x,
                                 long n,
                                 long n1,
                                 long n2)
    :w_(0)
    ,pw_(0)
    ,ifault_(0) 
{
    swilk(false, x, n, n1, n2, a_, w_, pw_, ifault_);
}

ms_shapiro_wilk::~ms_shapiro_wilk()
{
}


double ms_shapiro_wilk::getResult() const
{
    return w_;
}


double ms_shapiro_wilk::getPValue() const
{
    return pw_;
}


double ms_shapiro_wilk::getErrorCode() const
{
    return ifault_;
}


void ms_shapiro_wilk::appendSampleValue(double x) {
    std::deque<std::pair<size_t,double> >::size_type idx = x_.size();
    x_.push_back(std::make_pair(idx, x));
    w_      = 0;
    pw_     = 0;
    ifault_ = 0;
}

void ms_shapiro_wilk::clearSampleValues() {
    x_.clear();
    w_      = 0;
    pw_     = 0;
    ifault_ = 0;
}

void ms_shapiro_wilk::calculate(long n, long n1, long n2) {
    swilk(false, x_, n, n1, n2, a_, w_, pw_, ifault_);

}


void ms_shapiro_wilk::swilk(bool init,
                            std::deque<std::pair<size_t,double> > x,
                            long n,
                            long n1,
                            long n2,
                            std::deque<double> &a,
                            double &w,
                            double &pw,
                            int &ifault) {
    /* Algorithm AS R94, Journal of the Royal Statistical Society Series C
     * (Applied Statistics) vol. 44, no. 4, pp. 547-551 (1995).
     */

    w = 0; // stop the behaviour change 'feature'
    pw = 0;
    ifault = 0;

    /* Initialized data */
    static const double zero = 0.0;
    static const double one = 1.0;
    static const double two = 2.0;
    static const double three = 3.0;

    static const double z90 = 1.2816;
    static const double z95 = 1.6449;
    static const double z99 = 2.3263;
    static const double zm = 1.7509;
    static const double zss = 0.56268;
    static const double bf1 = 0.8378;
    static const double xx90 = 0.556;
    static const double xx95 = 0.622;
    static const double sqrth = 0.70711;/* sqrt(1/2) = .7071068 */
    static const double small_value = 1e-19;
    static const double pi6 = 1.909859;
    static const double stqr = 1.047198;

    /* polynomial coefficients */
    static const double g[2] = { -2.273,.459 };
    static const double c1[6] = { 0.,.221157,-.147981,-2.07119, 4.434685, -2.706056 },
                     c2[6] = { 0.,.042981,-.293762,-1.752461,5.682633, -3.582633 };
    static const double c3[4] = { .544,-.39978,.025054,-6.714e-4 };
    static const double c4[4] = { 1.3822,-.77857,.062767,-.0020322 };
    static const double c5[4] = { -1.5861,-.31082,-.083751,.0038915 };
    static const double c6[3] = { -.4803,-.082676,.0030302 };
    static const double c7[2] = { .164,.533 };
    static const double c8[2] = { .1736,.315 };
    static const double c9[2] = { .256,-.00635 };

    /* System generated locals */
    double r__1;

    a.resize(n, 0.0);

    /*
     * Auxiliary routines : poly()  {below}
     */

    /* Local variables */
    int i, j, ncens, i1, nn2;

    double zbar, ssassx, summ2, ssumm2, gamma, delta, range;
    double a1, a2, an, bf, ld, m, s, sa, xi, sx, xx, y, w1;
    double fac, asa, an25, ssa, z90f, sax, zfm, z95f, zsd, z99f, rsn, ssx, xsx;

    /* Parameter adjustments */
    pw = 1.0;
    ifault = 0;
    if (w >= 0.0) {
        w = 1.0;
    }
    an = (double) (n);
    nn2 = n / 2;
    if (n2 < nn2) {
        ifault = 3; 
        return;
    }
    if (n < 3) {
        ifault = 1; 
        return;
    }

    /*  If INIT is false, calculate coefficients a[] for the test */
    if (! (init)) {
        if (n==3) {
                a[0] = sqrth;
        } else {
            an25 = an + .25;
            summ2 = zero;
            for (i = 1; i <= n2; ++i) {
                a[i-1] = (double) ppnd7((i - .375f) / an25, ifault);
                // a[i-1] = (double) this->ppnd7((i - .375f) / an25, ifault);

                if (ifault) {
                    ifault = 8;
                    return;
                }
                summ2 += (a[i-1] * a[i-1]);
            }
            summ2 *= two;
            ssumm2 = std::sqrt(summ2);
            rsn = one / std::sqrt(an);
            a1 = poly(c1, 6, rsn) - a[0] / ssumm2;

            /* Normalize a[] */
            if (n > 5) {
                i1 = 3;
                a2 = -a[1] / ssumm2 + poly(c2, 6, rsn);
                fac = std::sqrt((summ2 - two * a[0] * a[0] - two * a[1] * a[1])
                            / (one - two * a1 * a1 - two * a2 * a2));
                a[1] = a2;
            } else {
                i1 = 2;
                fac = std::sqrt((summ2 - two * a[0] * a[0]) / ( one  - two * a1 * a1));
            }
            a[0] = a1;
            for (i = i1; i <= nn2; ++i) a[i-1] /= - fac;
        }
        init = true;
    }
    if (n1 < 3) {
        ifault = 1; 
        return;
    }
    ncens = n - n1;
    if (ncens < 0 || (ncens > 0 && n < 20)) {
        ifault = 4; 
        return;
    }
    delta = (double) ncens / an;
    if (delta > 0.8) {
        ifault = 5; 
        return;
    }

    /*  If W input as negative, calculate significance level of -W */
    if (w < zero) {
        w1 = 1. + w;
        ifault = 0;
        goto L70;
    }

    /*  Check for zero range */
    range = x[n1 - 1].second - x[0].second;
    if (range < small_value) {
        ifault = 6; 
        return;
    }

    /*  Check for correct sort order on range - scaled X */

    /* *ifault = 7; <-- a no-op, since it is set 0, below, in ANY CASE! */
    ifault = 0;
    xx = x[0].second / range;
    sx = xx;
    sa = -a[0];
    j = n - 1;
    for (i = 2; i <= n1; ++i) {
        xi = x[i-1].second / range;
        if (xx - xi > small_value) {
            /* Fortran had:  print *, "ANYTHING"
            * but do NOT; it *does* happen with sorted x (on Intel GNU/linux):
            *  shapiro.test(c(-1.7, -1,-1,-.73,-.61,-.5,-.24, .45,.62,.81,1))
            */
            ifault = 7;
            return;
        }
        sx += xi;
        if (i != j) sa += sign(1,i - j) * a[std::min(i,j)-1];
        xx = xi;
        --j;
    }
    if (n > 5000) {
        ifault = 2;
    }

/*  Calculate W statistic as squared correlation
    between data and coefficients */

    sa /= n1;
    sx /= n1;
    ssa = ssx = sax = zero;
    j = n;
    for (i = 1; i <= n1; ++i, --j) {
        if (i != j)
            asa = sign(1,i - j) * a[std::min(i,j)-1] - sa;
        else
            asa = -sa;
        xsx = x[i-1].second / range - sx;
        ssa += asa * asa;
        ssx += xsx * xsx;
        sax += asa * xsx;
    }

/*  W1 equals (1-W) claculated to avoid excessive rounding error
    for W very near 1 (a potential problem in very large samples) */

    ssassx = std::sqrt(ssa * ssx);
    w1 = (ssassx - sax) * (ssassx + sax) / (ssa * ssx);
L70:
    w = 1. - w1;

/*  Calculate significance level for W */

    if (n == 3) {/* exact P value : */
        pw = pi6 * (std::asin(std::sqrt(w)) - stqr);
        return;
    }
    y = std::log(w1);
    xx = std::log(an);
    m = zero;
    s = one;
    if (n <= 11) {
        gamma = poly(g, 2, an);
        if (y >= gamma) {
            pw = small_value;/* FIXME: rather use an even small_valueer value, or NA ? */
            return;
        }
        y = -std::log(gamma - y);
        m = poly(c3, 4, an);
        s = std::exp(poly(c4, 4, an));
    } else {/* n >= 12 */
        m = poly(c5, 4, xx);
        s = std::exp(poly(c6, 3, xx));
    }
    /*DBG printf("c(w1=%g, w=%g, y=%g, m=%g, s=%g)\n",w1,*w,y,m,s); */

    if (ncens > 0) {/* <==>  n > n1 */
    /*  Censoring by proportion NCENS/N.
        Calculate mean and sd of normal equivalent deviate of W. */

        ld = -std::log(delta);
        bf = one + xx * bf1;
        r__1 = pow(xx90, (double) xx);
        z90f = z90 + bf * std::pow(poly(c7, 2, r__1), (double) ld);
        r__1 = pow(xx95, (double) xx);
        z95f = z95 + bf * std::pow(poly(c8, 2, r__1), (double) ld);
        z99f = z99 + bf * std::pow(poly(c9, 2, xx), (double)ld);

        /* Regress Z90F,...,Z99F on normal deviates Z90,...,Z99 to get
         * pseudo-mean and pseudo-sd of z as the slope and intercept 
         */

        zfm = (z90f + z95f + z99f) / three;
        zsd = (z90 * (z90f - zfm) +
                    z95 * (z95f - zfm) + z99 * (z99f - zfm)) / zss;
        zbar = zfm - zsd * zm;
        m += zbar * s;
        s *= zsd;
    }
    pw = alnorm((y-m)/s, true);
    /*  = alnorm_(dble((Y - M)/S), 1); */

    // Results are returned in w, pw and ifault
    return;
} /* swilk */


double ms_shapiro_wilk::poly(const double *cc, int nord, double x)
{
    /* Auxiliary procedure in algorithm AS 181, Journal of the Royal 
     * Statistical Society Series C (Applied Statistics) vol. 31, no. 2, 
     * pp. 176-180 (1982).
     */

    /* Local variables */

    int n2,j;
    double ret_val = cc[0];
    if (nord == 1) return ret_val;
    double p = x * cc[nord-1];
    if (nord != 2) {
        n2 = nord - 2;
        j = n2 + 1;
        for (int i=1;i<=n2;i++) {
            p = (p + cc[j-1])*x;
            j--;
        }
    }
    ret_val = ret_val + p;
    return ret_val;
} /* poly */


double ms_shapiro_wilk::ppnd7(double p, int &ifault) {
    /* Algorithm AS 241, Journal of the Royal Statistical Society Series C
     * (Applied Statistics) vol. 26, no. 3, pp. 118-121 (1977).
     */
  
    static const double zero = 0.0;
    static const double one = 1.0;
    static const double half = 0.5;
    static const double split1 = 0.425;
    static const double split2 = 5.0;
    static const double const1 = 0.180625;
    static const double const2 = 1.6;
    static const double a0 = 3.3871327179E+00;
    static const double a1 = 5.0434271938E+01;
    static const double a2 = 1.5929113202E+02;
    static const double a3 = 5.9109374720E+01;
    static const double b1 = 1.7895169469E+01;
    static const double b2 = 7.8757757664E+01;
    static const double b3 = 6.7187563600E+01;
    static const double c0 = 1.4234372777E+00;
    static const double c1 = 2.7568153900E+00;
    static const double c2 = 1.3067284816E+00;
    static const double c3 = 1.7023821103E-01;
    static const double d1 = 7.3700164250E-01;
    static const double d2 = 1.2021132975E-01;
    static const double e0 = 6.6579051150E+00;
    static const double e1 = 3.0812263860E+00;
    static const double e2 = 4.2868294337E-01;
    static const double e3 = 1.7337203997E-02;
    static const double f1 = 2.4197894225E-01;
    static const double f2 = 1.2258202635E-02;

    double normal_dev;
    double q;
    double r;

    ifault = 0;
    q = p - half;
    if (std::abs(q) <= split1) {
        r = const1 - q * q;
        normal_dev = q * (((a3 * r + a2) * r + a1) * r + a0) / 
                     (((b3 * r + b2) * r + b1) * r + one);
        return normal_dev;
    } else {
        if (q < zero) {
            r = p;
        } else {
            r = one - p;
        }
        if (r <= zero) {
            ifault = 1;
            normal_dev = zero;
            return normal_dev;
        }
        r = std::sqrt(-std::log(r));
        if (r <= split2) {
            r = r - const2;
            normal_dev = (((c3 * r + c2) * r + c1) * r + c0) / ((d2 * r + d1) * r + one);
        } else {
            r = r - split2;
            normal_dev = (((e3 * r + e2) * r + e1) * r + e0) / ((f2 * r + f1) * r + one);
        }
        if (q < zero) { normal_dev = - normal_dev; }
        return normal_dev;
    }
}


double ms_shapiro_wilk::alnorm(double x,bool upper) {
    /* Algorithm AS 66, Journal of the Royal Statistical Society Series C
     * (Applied Statistics) vol. 22, pp. 424-427 (1973).
     */

    static const double zero = 0;
    static const double one = 1;
    static const double half = 0.5;
    static const double con = 1.28;
    static const double ltone = 7.0;
    static const double utzero = 18.66;
    static const double p = 0.398942280444;
    static const double q = 0.39990348504;
    static const double r = 0.398942280385;   
    static const double a1 = 5.75885480458;
    static const double a2 = 2.62433121679;
    static const double a3 = 5.92885724438;  
    static const double b1 = -29.8213557807;
    static const double b2 = 48.6959930692;
    static const double c1 = -3.8052E-8;
    static const double c2 = 3.98064794E-4;
    static const double c3 = -0.151679116635;
    static const double c4 = 4.8385912808;
    static const double c5 = 0.742380924027;
    static const double c6 = 3.99019417011;  
    static const double d1 = 1.00000615302;
    static const double d2 = 1.98615381364;
    static const double d3 = 5.29330324926;  
    static const double d4 = -15.1508972451;
    static const double d5 = 30.789933034;
  
    double alnorm;
    double z;
    double y;
    bool up = upper;
    z = x;
    if (z < zero){
        up = !up;
        z = -z;
    }
    if (z <= ltone || (up && z <= utzero)) {
        y = half * z * z;
        if (z > con) {
            alnorm = r * std::exp(-y) / (z + c1 + d1 / (z + c2 + d2 / (z + c3 + d3 
                    / (z + c4 + d4 / (z + c5 + d5 / (z + c6))))));
        } else {
            alnorm = half - z * (p - q * y/(y + a1 + b1 /(y + a2 + b2/(y + a3))));
        }
    } else {
        alnorm = zero;
    }

    if(!up) { alnorm = one - alnorm; }
    return alnorm;
}

long ms_shapiro_wilk::sign(long x,long y) {
    if (y<0) 
        return -std::labs(x);
    else
        return std::labs(x);
}
void get_n1_n2(std::vector<double>  &data,long*n1,long*n2){

    // Calculate the mean and variance of the sorted data
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    double variance = std::accumulate(data.begin(), data.end(), 0.0, [mean](double acc, double x) {
        return acc + (x - mean) * (x - mean);
    }) / (data.size() - 1);

    // Calculate the differences between each sorted data point and the mean
    std::vector<double> diff(data.size());
    std::transform(data.begin(), data.end(), diff.begin(), [mean](double x) {
        return x - mean;
    });

    // Split the sorted data into two equal-sized groups
    std::vector<double> group1(diff.begin(), diff.begin() + data.size() / 2);
    std::vector<double> group2(diff.begin() + data.size() / 2, diff.end());

    // Calculate n1 and n2
    double sum1 = std::accumulate(group1.begin(), group1.end(), 0.0, [](double acc, double x) {
        return acc + x * x;
    });
    double sum2 = std::accumulate(group2.begin(), group2.end(), 0.0, [](double acc, double x) {
        return acc + x * x;
    });
    *n1 = static_cast<long>(sum1 / variance + 0.5);
    *n2 = static_cast<long>(sum2 / variance + 0.5);

}
int main() {
    std::vector<double> data {
        1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,
        1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,
        1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0,1.0, 2.0, 3.0, 4.0, 5.0
    };
    std::sort(data.begin(), data.end());

    std::deque<std::pair<size_t, double>> x;
    for (size_t i = 0; i < data.size(); i++) {
        x.push_back(std::make_pair(i+1, data[i]));
    }
    int n = x.size();
    long n1=0 ;
    long n2 =0;
    get_n1_n2(data,&n1,&n2);
    
    std::cout << "n : " << n << std::endl;
    std::cout << "n1 : " << n1 << std::endl;
    std::cout << "n2 :" << n2 << std::endl;
    // Create ms_shapiro_wilk object and calculate W and p-value
    ms_shapiro_wilk sw(x, x.size(), n1, n2);

    // sw.calculate(x.size(), 0, 0);
    // sw.swilk(true, x, x.size(), 0, 0, sw.a_, sw.w_, sw.pw_, sw.ifault_);

    double W = sw.getResult();
    double p_value = sw.getPValue();
    double e = sw.getErrorCode();
    std::cout << "W statistic: " << W << std::endl;
    std::cout << "p-value: " << p_value << std::endl;
    std::cout << "getErrorCode :" << e << std::endl;

    return 0;

}

