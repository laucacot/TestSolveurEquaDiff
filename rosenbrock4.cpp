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


value_type k6 (value_type Te) //K6 SiH4 + e -> SiH3 + H + e
{
    value_type K6;
    K6=1.83E-9*pow((Te),-1)*exp(-(10.68)/(Te));
    return K6;

}

struct nsystem
{
    void operator()( const state_type &n , state_type &dndt , const value_type &t ) 
    { /*0=SiH3  1=SiH4*/
        dndt( 0 ) = k6(Te)*n_e*n(1);
        dndt( 1 ) =-k6(Te)*n_e*n(1) ;
    }
 value_type Te;
};

struct jacobian
{
    void operator()( const state_type &n , matrix_type &jacobi , const value_type &t , state_type &dfdt ) const
    {
        jacobi( 0 , 0 ) = 0;
        jacobi( 0 , 1 ) = k6(Te)*n_e;
        jacobi( 1 , 0 ) = 0;
        jacobi( 1 , 1 ) = -k6(Te)*n_e;
        dfdt( 0 ) = 0.0;
        dfdt( 1 ) = 0.0;
    }
value_type Te;
};

struct etemperature
{
  value_type operator()(value_type const& Te)
  {
    return -DP/n_e
    
    +k6(Te)*n[1]*10.68;
  }

  state_type n;
};


void write_density( const value_type t, const value_type Te, const state_type &n)
{
    cout << t << '\t' <<Te <<'\t' 
              << n[0] << '\t' << n[1] << '\t'
              << endl;
}

int main(int argc, char **argv)
{
cout <<"t"<<'\t'<<"Te"<<'\t'<<"SiH3"<<'\t'<<"SiH4"<<endl;

// Time variables
  value_type t = 0.0;
  value_type dt = 1.0e-8;
  value_type Tmax = 20.e-6;
  value_type NT = Tmax/dt;

  // Root finding variables
  value_type min = 0.001;
  value_type max = 20.0;
  boost::uintmax_t max_iter = 500;
 eps_tolerance<value_type> tol(30);

  // initial values
  value_type Te = 3.0;

  // Density vectors and initial condition
  state_type n_ini(Nbr_espece, 0.0); // initial conditions
  n_ini[0] = 0;
  n_ini[1] = n_SiH4_ini;
 
  state_type n_new(Nbr_espece, 0.0);  // first step same as initial conditions
  n_new = n_ini;
  state_type n_err(Nbr_espece, 0.0); //error

  // declare the functor etemperature
  etemperature etemp;
  // assign initial values to functor etemp
  etemp.n = n_ini;

  //cerr << "\n[ii] Electrons  = " << etemp.n[0] << endl;
  //cerr << "\n[ii] Metastables  = " << etemp.n[1] << endl;
//   cout << "\n[ii] n  = " << etemp.n[2] << endl;
 
//t1=clock();
  // Find Te first calculation
 // pair<value_type, value_type> pair_Te =\
                toms748_solve(etemp, min, max, tol, max_iter);

  //Te = pair_Te.first;
  cerr << "\n[ii] Initial Temperature  = " << Te << endl;
//t2=clock()-t1;
//cout<<"timesec"<<(float )t2/CLOCKS_PER_SEC << endl;

  // declare system and jacobian
  nsystem sys;
  jacobian jac;

  // declare stepper Rosenbrock
  stepper_type stepper;

  for (int i = 0; i <= NT; i++)
  {
    // update Te in system and jacobian
    sys.Te = Te;
    jac.Te = Te;

    // Integrate at least one step dt
    stepper.do_step( std::make_pair( sys, jac ), n_new, t, dt, n_err);

    // assign values to functor etemp
    etemp.n = n_new;
if (i%((int)(NT/100))==0)
{
    write_density(t, Te, n_new);
}
    // Find new Te
    //pair<value_type, value_type> pair_Te =\
                  toms748_solve(etemp, min, max, tol, max_iter);

    //Te = pair_Te.first;

    t+= dt;
    n_ini = n_new;//update
  }

/*float charge= (n_new[20]+n_new[4]+n_new[10]-n_new[0]-n_new[2]-n_new[3]-n_new[13]-n_new[15]-n_new[16]-n_new[17])/n_Arp_ini;

cerr<<"charge/dArp="<<charge<<endl;*/

float Si=(n_new[0]+n_new[1])/n_SiH4_ini;

cerr<<"Si="<<Si<<endl;

float H=(3*n_new[0]+4*n_new[1])/(4*n_SiH4_ini);

cerr<<"H="<<H<<endl;
  return 0;

}  
