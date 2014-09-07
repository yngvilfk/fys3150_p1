
#include <istream>
#include <iostream>
#include <print.h>




Matlab::Print::Print(const arma::Col<double>& vec1,
                     const arma::Col<double>& vec2,
                     const std::string&       filename)
{
    std::ofstream writeMatlab;
    writeMatlab.open(filename.c_str());

    // Write vec1 to matlab vector
    writeMatlab << "vector_1 = [" << "\n";

    for (uint i = 0; i < vec1.size(); ++i)
    {
        writeMatlab << vec1(i) << "\n";
    }
    writeMatlab << "];" << "\n";

    // Write vec2 to matlab vector
    writeMatlab << "vector_2 = [" << "\n";

    for (uint i = 0; i < vec2.size(); ++i)
    {
        writeMatlab << vec2(i) << "\n";
    }
    writeMatlab << "];" << "\n";

    writeMatlab.close();
}



Matlab::Print::~Print()
{
}
