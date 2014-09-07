#include <print.h>
#include <string>
#include <iostream> //i standardbiblioteket
#include <vector> //ligger i standardbiblioteket
#include <cstdlib>
#include <cmath>
#include <armadillo>
#include <time.h>
#include <fstream>

using namespace arma;
using namespace std;



namespace
{
   namespace Local
   {
      bool input(int,
                 char**,
                 int&         n,
                 std::string& fileName);

      void writeMatlab(std::ofstream& write,
                       arma::Col<double>&,
                       const std::string& name);
   }
}



int main(int argv, char* argc[] )
{
   int n;
   std::string fileName;
   const bool runLU = Local::input(argv, argc, n, fileName);

   std::ofstream write(fileName.c_str());

   static_cast<void>(argv);



   //declarattion of variables
   const double h = 1.0/(static_cast<double>(n+1));
   arma::Col<double> f(n+2),
         v(n+2),       //numerical approximation
         temp(n+1),
         x(n+2),       //x-values
         u(n+2) ;      //analytical solution
   double j;

   clock_t start, finish;



   for(int i=0; i<=n+1 ; ++i)
   {
      j = static_cast<double>(i);
      f(i) = std::pow(h, 2)*100*std::exp(-10*j*h);         //f = h^2*f(x) = h^2*100exp(10x)
      u(i) = 1-(1- std::exp(-10))*j*h - std::exp(-10*j*h); // analytical solution
      x(i) = j*h;   //x= 0, h, 2h,...,jh,..., nh, 1
   }

   //solve  -v''(x) = f(x)
   double b = 2.0;
   double c=-1.0;

   temp.fill(1.0);          //all elements in temp is set to 1.0

   double btemp = b;
   v(0) = 0.0;
   v(n+1) = 0.0;
   v(1) = f(1)/btemp;

   start = clock();         //start timer

   //Decomposition and forward substitution
   for(int i=2 ; i <= n ; ++i)
   {
      temp(i) = c/btemp;
      btemp = b + temp(i);
      v(i) = (f(i) + v(i-1))/btemp;
   }
   //backsubstitution
   for(int i = n-1 ; i >= 1 ; --i)
   {
      v(i) -= temp(i+1)*v(i+1);
      //        std::cout << "analytical u: " << u(i) << " numerical u: " << v(i) << std::endl;
   }
   finish = clock();        //stop timer

   cout << "time tridiagonal algorithm: " << static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC ) << " s" << endl;

   //write to file for plotting
   const std::string fileName_v("v_vec");
   Local::writeMatlab(write, v, fileName_v);

   const std::string fileName_u("u_vec");
   Local::writeMatlab(write, u, fileName_u);

   const std::string fileName_x("x_vec");
   Local::writeMatlab(write, x, fileName_x);


   //create matrix A
   if (runLU) //only runs if n<=10000
   {
      mat A(n,n);
      A.fill(0);
      A(0,0)=2.0;
      A(n-1,n-1)=2.0;
      A(0,1)= -1.0;
      A(n-1,n-2) = -1;

      for(int i=1; i<n-1 ;++i)
      {
         A(i,i) = 2.0;
         A(i,i+1) = -1;
         A(i,i-1) = -1;
      }

      vec F
            = f.subvec(1,n);        //make a new vector with right dimension

      start = clock();   //start timer



      //LU-decomposition, solve Av=LUv=f
      mat L, U, P;
      lu(L, U, P, A);

        //testing the LU-decomposition
       mat B = P.t()*L*U;

      //    cout << B << endl;

      //solve Ly=f and Uu=y to find u:
      vec y = solve(L, F);
      vec uLU = solve(U, y);
      finish = clock();       //stop timer

      cout << "time LU-decomposition: " << static_cast<double>(finish - start)/static_cast<double>(CLOCKS_PER_SEC) << " s" << endl;
    }


   // compute relative error
   arma::Col<double> epsilon(n+1);
   for (int i = 1; i <= n; ++i)
   {
      epsilon(i) = log10(fabs((v(i)-u(i))/u(i)));
   }
   cout << "log10(epsilon) for n = " << n << ": " << arma::min(epsilon) << endl;
   //reative error should be constant. This can be proved


   write.close();
   return 0;
}


//functions to write to read from terminal and write to file
bool
Local::input(int         argv,
             char**      argc,
             int&         n,
             std::string& fileName)
{
   if (argv == 3)
   {
      n = atoi(argc[1]);
      fileName.append(argc[2]);
   }
   else if (argv == 2)
   {
      std::cout << "filename: ";
      std::cin >> fileName;
      std::cout << std::endl;
      n = atoi(argc[1]);
   }
   else if (argv <= 1)
   {
      std::cout << "n       : ";
      std::cin >> n;
      std::cout << std::endl;

      std::cout << "filename: ";
      std::cin >> fileName;
      std::cout << std::endl;
   }

   if (n <= 1000)
      return true;
   else
      return false;
}


void
Local::writeMatlab(std::ofstream& write,
                   arma::Col<double>& vec,
                   const std::string& name)
{
   write << name << " = [" << "\n";
   for (std::size_t iVec = 0; iVec < vec.size(); ++iVec)
   {
      write << vec(iVec) << "\n";
   }

   write << "]" << "\n\n";
}
