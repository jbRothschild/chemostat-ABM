// sampled_distribution.hh
#ifndef SAMPLED_DISTRIBUTION
#define SAMPLED_DISTRIBUTION

#include <algorithm>
#include <vector>
#include <random>
#include <stdexcept>

template <typename T = double, bool Interpolate = true>
class Sampled_distribution
{
public:
    using CDFFunc = T (*)(T);

    Sampled_distribution(CDFFunc cdfFunc, T low, T high, unsigned resolution = 200)
        : mLow(low), mHigh(high), mRes(resolution), mDist(0.0, 1.0)
    {
        if (mLow >= mHigh) throw InvalidBounds();

        mSampledCDF.resize(mRes + 1);
        const T cdfLow = cdfFunc(low);
        const T cdfHigh = cdfFunc(high);
        T last_p = 0;
        for (unsigned i = 0; i < mSampledCDF.size(); ++i) {
            const T x = i/mRes*(mHigh - mLow) + mLow;
            const T p = (cdfFunc(x) - cdfLow)/(cdfHigh - cdfLow); // normalising
            if (! (p >= last_p)) throw CDFNotMonotonic();
            mSampledCDF[i] = Sample{p, x};
            last_p = p;
        }
    }

    template <typename Generator>
    T operator()(Generator& g)
    {
        T cdf = mDist(g);
        auto s = std::upper_bound(mSampledCDF.begin(), mSampledCDF.end(), cdf);
        auto bs = s - 1;
        if (Interpolate && bs >= mSampledCDF.begin()) {
            const T r = (cdf - bs->prob)/(s->prob - bs->prob);
            return r*bs->value + (1 - r)*s->value;
        }
        return s->value;
    }

private:
    struct InvalidBounds : public std::runtime_error { InvalidBounds() : std::runtime_error("") {} };
    struct CDFNotMonotonic : public std::runtime_error { CDFNotMonotonic() : std::runtime_error("") {} };

    const T mLow, mHigh;
    const double mRes;

    struct Sample {
        T prob, value;
        friend bool operator<(T p, const Sample& s) { return p < s.prob; }
    };

    std::vector<Sample> mSampledCDF;
    std::uniform_real_distribution<> mDist;
};

#endif
