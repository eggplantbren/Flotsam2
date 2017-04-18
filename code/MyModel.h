#ifndef Flotsam2_MyModel
#define Flotsam2_MyModel

#include "DNest4/code/DNest4.h"
#include <ostream>
#include "celerite/celerite.h"

namespace Flotsam2
{

class MyModel
{
    private:

        // Hyperparameters for mean levels of images
        double mu_magnitudes;
        double sig_magnitudes;

        // Normals for mean levels of images
        std::vector<double> magnitude_ns;
        std::vector<double> magnitudes;
        void compute_magnitudes();

        // Amplitude and timescale for QSO variability
        double qso_amplitude;
        double qso_timescale;

        // Time delays (relative to image 0)
        std::vector<double> time_delays;

        // Amplitudes and timescales for extra
        // correlated residual ("microlensing")
        double mu_amplitude, sig_log_amplitude;
        double mu_timescale, sig_log_timescale;
        std::vector<double> amplitude_ns;
        std::vector<double> timescale_ns;

        // Error bar boost parameter
        double u_boost;
        double sigma_boost_factor;
        void compute_sigma_boost_factor();

    public:
        // Constructor only gives size of params
        MyModel();

        // Generate the point from the prior
        void from_prior(DNest4::RNG& rng);

        // Metropolis-Hastings proposals
        double perturb(DNest4::RNG& rng);

        // Likelihood function
        double log_likelihood() const;

        // Print to stream
        void print(std::ostream& out) const;

        // Return string with column information
        std::string description() const;
};

} // namespace Flotsam2

#endif

