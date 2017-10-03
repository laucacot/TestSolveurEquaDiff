#include <iostream>
#include <string>
#include <cmath>
#include <vector>
// Optimization
#define BOOST_UBLAS_NDEBUG
#include <fstream>
#include <exception>
#include <string>
#include <utility>
#include <boost/array.hpp>
#include <boost/numeric/odeint.hpp>
#include <time.h>
#include <stdlib.h>








using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::math::tools;

// type definitions
typedef double value_type;// or typedef float value_type;
typedef boost::numeric::ublas::vector< value_type > state_type;
typedef boost::numeric::ublas::matrix< value_type > matrix_type;
typedef rosenbrock4< value_type > stepper_type;

// constants
const value_type pressure = 13.332237; //pascal soit 0.1 torr
const value_type Tg = 0.02758;//320 K
const value_type L = 3e-2; //distance netre deux plaques en m
const value_type k_b = 1.38064852e-23; // constante de boltzman en J/K
const value_type K = 8.6173303e-5; //constqnte de boltzman en eV/k
const value_type Da_Arp=0.005; //diffusion des Ar+ en m2/s
const value_type pi = M_PI;
const value_type diff = pow((pi/L), 2);
const value_type DP = 2e19;//puissance totale du systeme
const value_type n_Ar = (0.1/760)*2.69e25;//pressure/(320*k_b);
const value_type n_SiH4_ini = n_Ar/100.;
const value_type n_e = 1.e16;
const int Nbr_espece=2;




  const   float k6=1.83E-9*pow((3),-1)*exp(-(10.68)/(3));



    void systeme ( const state_type &n , state_type &dndt , double t ) 
    { /*0=SiH3  1=SiH4*/
        dndt( 0 ) = k6*n_e*n(1);
        dndt( 1 ) =-k6*n_e*n(1) ;
    }





int write_density (const  state_type &n , const double t )
{
    
cout << t <<'\t' 
              << n[0] << '\t' << n[1] 
              << endl;
return 0;
}

int main(int argc, char **argv)
{
cout <<"t"<<'\t'<<"Te"<<'\t'<<"SiH3"<<'\t'<<"SiH4"<<endl;

// Time variables
  double t = 0.0;
  double dt = 1.0e-8;
  double Tmax = 20.e-6;
  value_type NT = Tmax/dt;

 
  // initial values
  value_type Te = 3.0;

  // Density vectors and initial condition
  state_type n_ini={0,n_SiH4_ini}; // initial conditions
 
  
  cerr << "\n[ii] Initial Temperature  = " << Te << endl;

  
    // Integrate at least one step dt
   integrate (systeme,n_ini, 0, 5,0.1,write_density);

  



float Si=(n_ini[0]+n_ini[1])/n_SiH4_ini;

cerr<<"Si="<<Si<<endl;

float H=(3*n_ini[0]+4*n_ini[1])/(4*n_SiH4_ini);

cerr<<"H="<<H<<endl;
  return 0;

}  
