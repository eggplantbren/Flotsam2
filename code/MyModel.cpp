#include "Data.h"
#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <sstream>

namespace Flotsam2
{

MyModel::MyModel()
:magnitude_ns(Data::get_instance().get_num_images())
,magnitudes(Data::get_instance().get_num_images())
{

}

void MyModel::from_prior(DNest4::RNG& rng)
{
    // Magnitudes could be relative or not, but this should
    // cover both cases pretty well
    mu_magnitudes = 30.0*rng.randn();
    sig_magnitudes = exp(log(0.1) + log(100.0)*rng.rand());

    for(size_t i=0; i<magnitude_ns.size(); ++i)
    {
        magnitude_ns[i] = rng.randn();
        magnitudes[i] = mu_magnitudes + sig_magnitudes*magnitude_ns[i];
    }
}

void MyModel::compute_sigma_boost_factor()
{
    if(u_boost <= 0.5)
    {
        sigma_boost_factor = 1.0;
    }
    else
    {
        double u = 2*(u_boost - 0.5);
        sigma_boost_factor = exp(-log(1.0 - u));
    }
}

double MyModel::perturb(DNest4::RNG& rng)
{
    double logH = 0.0;

    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    return logL;
}

void MyModel::print(std::ostream& out) const
{

}

std::string MyModel::description() const
{
    std::stringstream s;

    return s.str();
}

} // namespace Flotsam2

