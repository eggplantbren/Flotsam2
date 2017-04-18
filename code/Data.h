#ifndef Flotsam2_Data
#define Flotsam2_Data

#include <vector>
#include "Eigen/Dense"

namespace Flotsam2
{

class Data
{
    private:
        std::vector<double> t, y, sig;
        std::vector<size_t> image;
        Eigen::VectorXd tt, yy, var;
        size_t num_images;

    public:
        Data();
        void load(const char* filename);

        // Getters
        const std::vector<double>& get_t()     const { return t; }
        const std::vector<double>& get_y()     const { return y; }
        const std::vector<double>& get_sig()   const { return sig; }
        const std::vector<size_t>& get_image() const { return image; }
        size_t get_num_images() const { return num_images; }

        // Getters of eigen vectors
        const Eigen::VectorXd& get_tt()  const { return tt; }
        const Eigen::VectorXd& get_yy()  const { return yy; }
        const Eigen::VectorXd& get_var() const { return var; }

        // Summaries
        double get_t_range() const { return (t.back() - t[0]); }

    // Singleton
    private:
        static Data instance;
    public:
        static Data& get_instance() { return instance; }
};

} // namespace

#endif

