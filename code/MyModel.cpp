#include "Data.h"
#include "MyModel.h"
#include "DNest4/code/DNest4.h"
#include <sstream>

namespace Flotsam2
{

MyModel::MyModel()
:magnitude_ns(Data::get_instance().get_num_images())
,magnitudes(Data::get_instance().get_num_images())
,time_delays(Data::get_instance().get_num_images(), 0.0)
,amplitude_ns(Data::get_instance().get_num_images())
,timescale_ns(Data::get_instance().get_num_images())
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

    // Grab timescale from data
    double t_range = Data::get_instance().get_t_range();
    qso_amplitude = -log(1.0 - rng.rand());
    qso_timescale = exp(log(0.1*t_range) + log(100.0)*rng.rand());

    // Time delay prior
    double width = 0.1*t_range;
    DNest4::Cauchy cauchy(0.0, width);
    for(size_t i=1; i<time_delays.size(); ++i)
    {
        do
        {
            time_delays[i] = cauchy.generate(rng);
        }while(std::abs(time_delays[i]) > t_range);
    }

    // Microlensing-related hyperparameters
    mu_amplitude = -log(1.0 - rng.rand());
    sig_log_amplitude = 3.0*rng.rand();
    mu_timescale = exp(log(0.1*t_range) + log(100.0)*rng.rand());
    sig_log_timescale = 3.0*rng.rand();

    // Microlensing-related normals
    for(double& n: amplitude_ns)
        n = rng.randn();
    for(double& n: timescale_ns)
        n = rng.randn();

    u_boost = rng.rand();
    compute_sigma_boost_factor();
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

    // Grab timescale from data
    double t_range = Data::get_instance().get_t_range();

    bool bulk = rng.rand() <= 0.5;
    if(bulk)
    {
        int which = rng.rand_int(4);

        if(which == 0)
        {
            int i = rng.rand_int(magnitude_ns.size());

            logH -= -0.5*pow(magnitude_ns[i], 2);
            magnitude_ns[i] += rng.randh();
            logH += -0.5*pow(magnitude_ns[i], 2);

            compute_magnitudes();
        }
        else if(which == 1)
        {
            // Which time delay to change
            int i = 1 + rng.rand_int(time_delays.size() - 1);

            double width = 0.1*t_range;
            DNest4::Cauchy cauchy(0.0, width);
            logH += cauchy.perturb(time_delays[i], rng);
            if(std::abs(time_delays[i]) > t_range)
                return -1E300;
        }
        else if(which == 2)
        {
            int i = rng.rand_int(amplitude_ns.size());

            logH -= -0.5*pow(amplitude_ns[i], 2);
            amplitude_ns[i] += rng.randh();
            logH += -0.5*pow(amplitude_ns[i], 2);
        }
        else
        {
            int i = rng.rand_int(timescale_ns.size());

            logH -= -0.5*pow(timescale_ns[i], 2);
            timescale_ns[i] += rng.randh();
            logH += -0.5*pow(timescale_ns[i], 2);
        }
    }
    else
    {
        int which = rng.rand_int(9);

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
        else if(which == 2)
        {
            qso_amplitude = 1.0 - exp(-qso_amplitude);
            qso_amplitude += rng.randh();
            DNest4::wrap(qso_amplitude, 0.0, 1.0);
            qso_amplitude = -log(1.0 - qso_amplitude);
        }
        else if(which == 3)
        {
            // Grab timescale from data
            double t_range = Data::get_instance().get_t_range();
            qso_timescale = log(qso_timescale);
            qso_timescale += log(100.0)*rng.randh();
            DNest4::wrap(qso_timescale, log(0.1*t_range), log(10.0*t_range));
            qso_timescale = exp(qso_timescale);
        }
        else if(which == 4)
        {
            mu_amplitude = 1.0 - exp(-mu_amplitude);
            mu_amplitude += rng.randh();
            DNest4::wrap(mu_amplitude, 0.0, 1.0);
            mu_amplitude = -log(1.0 - mu_amplitude);
        }
        else if(which == 5)
        {
            sig_log_amplitude += 3.0*rng.randh();
            DNest4::wrap(sig_log_amplitude, 0.0, 3.0);
        }
        else if(which == 6)
        {
            mu_timescale = log(mu_timescale);
            mu_timescale += log(100.0)*rng.randh();
            DNest4::wrap(mu_timescale, log(0.1*t_range), log(10.0*t_range));
            mu_timescale = exp(mu_timescale);       
        }
        else if(which == 7)
        {
            sig_log_timescale += 3.0*rng.randh();
            DNest4::wrap(sig_log_timescale, 0.0, 3.0);
        }
        else
        {
            u_boost += rng.randh();
            DNest4::wrap(u_boost, 0.0, 1.0);
            compute_sigma_boost_factor();
        }
    }

    return logH;
}

double MyModel::log_likelihood() const
{
    double logL = 0.0;

    // Grab the data
    const Data& data = Data::get_instance();

    // Copies of, or references to, things in the data
    auto image = data.get_image();
    std::vector<double> t = data.get_t();
    std::vector<double> y = data.get_y();
    std::vector<double> sig = data.get_sig();

    // Adjust for time delays
    for(size_t i=0; i<t.size(); ++i)
        t[i] -= time_delays[image[i]];

    // Subtract the magnitudes
    for(size_t i=0; i<y.size(); ++i)
        y[i] -= magnitudes[image[i]];

    // Argsort the times (after adjustment for time delays)
    std::vector<size_t> indices = DNest4::argsort(t);

    // Create eigen vector versions of the data, with the
    // new ordering.
    Eigen::VectorXd tt(t.size());
    Eigen::VectorXd yy(t.size());
    Eigen::VectorXd var(t.size());
    double fsq = pow(sigma_boost_factor, 2);
    for(size_t i=0; i<t.size(); ++i)
    {
        tt(i) = t[indices[i]];
        yy(i) = y[indices[i]];
        var(i) = fsq*pow(sig[indices[i]], 2);
    }

    // QSO term
    Eigen::VectorXd alpha_real(1 + data.get_num_images()),
                    beta_real (1 + data.get_num_images());
    alpha_real(0) = qso_amplitude;
    beta_real(0) = 1.0/qso_timescale;

    // Microlensing terms
    for(size_t i=0; i<data.get_num_images(); ++i)
    {
        alpha_real(i+1) = mu_amplitude*
                            exp(sig_log_amplitude*amplitude_ns[i]);
        beta_real(i+1)  = 1.0/mu_timescale/
                            exp(sig_log_timescale*timescale_ns[i]);
    }

    // Celerite solver
    celerite::solver::BandSolver<double> solver(true);
    solver.compute(alpha_real, beta_real,
                   Eigen::VectorXd(0),
                   Eigen::VectorXd(0),
                   Eigen::VectorXd(0),
                   Eigen::VectorXd(0),
                   tt, var);

    logL += -0.5*log(2*M_PI)*yy.size();
    logL += -0.5*solver.log_determinant();
    logL += -0.5*solver.dot_solve(yy);

    return logL;
}

void MyModel::print(std::ostream& out) const
{
    out << mu_magnitudes << ' ';
    out << sig_magnitudes << ' ';

    for(double m: magnitudes)
        out << m << ' ';

    out << qso_amplitude << ' ';
    out << qso_timescale << ' ';

    for(double tau: time_delays)
        out << tau << ' ';

    out << mu_amplitude << ' ';
    out << sig_log_amplitude << ' ';
    out << mu_timescale << ' ';
    out << sig_log_timescale << ' ';

    for(double n: amplitude_ns)
        out << (mu_amplitude*exp(sig_log_amplitude*n)) << ' ';
    for(double n: timescale_ns)
        out << (mu_timescale*exp(sig_log_timescale*n)) << ' ';

    out << sigma_boost_factor << ' ';
}

std::string MyModel::description() const
{
    std::stringstream s;

    s << "mu_magnitudes, ";
    s << "sig_magnitudes, ";
    for(size_t i=0; i<magnitudes.size(); ++i)
        s << "magnitudes[" << i << "], ";

    s << "qso_amplitude, qso_timescale, ";

    for(size_t i=0; i<magnitudes.size(); ++i)
        s << "time_delays[" << i << "], ";

    s << "mu_amplitude, ";
    s << "sig_log_amplitude, ";
    s << "mu_timescale, ";
    s << "sig_log_timescale, ";

    for(size_t i=0; i<amplitude_ns.size(); ++i)
        s << "ml_amplitudes[" << i << "], ";
    for(size_t i=0; i<amplitude_ns.size(); ++i)
        s << "ml_timescales[" << i << "], ";

    s << "sigma_boost_factor, ";

    return s.str();
}

} // namespace Flotsam2

