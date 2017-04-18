#include "Data.h"
#include <exception>
#include <fstream>
#include <iostream>
#include <set>

namespace Flotsam2
{

Data Data::instance;

Data::Data()
{

}

void Data::load(const char* filename)
{
    std::fstream fin(filename, std::ios::in);
    if(!fin)
    {
        std::cerr<<"# Error. Couldn't open file "<<filename<<"."<<std::endl;
        return;
    }

    // Empty the vectors
    t.clear();
    y.clear();
    sig.clear();
    image.clear();

    double temp1, temp2, temp3;
    double temp4;
    while(fin>>temp1 && fin>>temp2 && fin>>temp3 && fin>>temp4)
    {
        t.push_back(temp1);
        y.push_back(temp2);
        sig.push_back(temp3);
        image.push_back((size_t)temp4);
    }
    std::cout<<"# Loaded "<<t.size()<<" data points from file "
            <<filename<<"."<<std::endl;
    fin.close();

    // Copy into eigen vectors
    tt.resize(t.size());
    yy.resize(y.size());
    var.resize(sig.size());
    for(size_t i=0; i<y.size(); ++i)
    {
        tt(i) = t[i];
        yy(i) = y[i];
        var(i) = sig[i] * sig[i];
    }

    // Use a set to count number of images
    std::set<size_t> temp;
    for(auto img: image)
        temp.insert(img);
    num_images = temp.size();

    std::cout<<"# Found data from "<<num_images<<" images."<<std::endl;
}

} // namespace Flotsam2

