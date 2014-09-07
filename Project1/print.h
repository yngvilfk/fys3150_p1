#ifndef PRINT_H
#define PRINT_H

#include <armadillo>
#include <fstream>
#include <memory>
#include <string>

namespace Matlab
{
   class Print
   {
   public:
      Print(const arma::Col<double>& vec1,
            const arma::Col<double>& vec2,
            const std::string&       filename);
      ~Print();

   protected:
      void writeFile(std::fstream&,
                     const arma::Col<double>&);

   private:
      typedef std::auto_ptr<arma::Col<double> > Vector;
      Vector vec_;
   };
}



#endif // PRINT_H
