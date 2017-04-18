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

    for(size_t i=0; i<magnitudes.size(); ++i)
        magnitude_ns[i] = rng.randn();

    compute_magnitudes();
}

void MyModel::compute_magnitudes()
{
    for(size_t i=0; i<magnitudes.size(); ++i)
        magnitudes[i] = mu_magnitudes + sig_magnitudes*magnitude_ns[i];
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

    int which = rng.rand_int(3);

    if(which == 0)
    {
        logH -= -0.5*pow(mu_magnitudes/30.0, 2);
        mu_magnitudes += 30.0*rng.randh();
        logH += -0.5*pow(mu_magnitudes/30.0, 2);

        compute_magnitudes();
    }
    else if(which == 1)
    {
        sig_magnitudes = log(sig_magnitudes);
        sig_magnitudes += log(100.0)*rng.randh();
        DNest4::wrap(sig_magnitudes, log(0.1), log(10.0));
        sig_magnitudes = exp(sig_magnitudes);

        compute_magnitudes();
    }
    else
    {
        int i = rng.rand_int(magnitude_ns.size());

        logH -= -0.5*pow(magnitude_ns[i], 2);
        magnitude_ns[i] += rng.randh();
        logH += -0.5*pow(magnitude_ns[i], 2);

        compute_magnitudes();
    }

    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    out << mu_magnitudes << ' ';
    out << sig_magnitudes << ' ';

    for(double m: magnitudes)
        out << m << ' ';
}

std::string MyModel::description() const
{
    std::stringstream s;

    s << "mu_magnitudes, ";
    s << "sig_magnitudes, ";
    for(size_t i=0; i<magnitudes.size(); ++i)
        s << "magnitudes[" << i << "], ";

    return s.str();
}

} // namespace Flotsam2

