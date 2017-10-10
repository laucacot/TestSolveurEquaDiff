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

#include <boost/phoenix/core.hpp>
#include <boost/phoenix/operator.hpp>

using namespace std;
using namespace boost::numeric::odeint;
using namespace boost::math::tools;
namespace phoenix = boost::phoenix;

// type definitions
typedef double value_type;// or typedef float value_type;
typedef vector< value_type > state_type;

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
const int Nbr_espece=3;
const int Te=3;

const int Nbr_K=1;
int jmax=Nbr_K;
int imax=9;
double **Tab;

value_type p1,p2,g1,g2,g3,g4,Tp,Tx,Tj;
state_type Kt(jmax, 0.0);
const float   k6=1.83E-9*pow((Te),-1)*exp(-(10.68)/(Te));


struct stiff_system
{
    inline
    void operator()( const state_type &n , state_type &dndt , const value_type &t ) const
    { /*0=SiH3  1=SiH4  2=H
        dndt( 0 ) = k6*n_e*n(1);
        dndt( 1 ) =-k6*n_e*n(1) ;
	dndt( 2 ) = k6*n_e*n(1);*/
    
for (int k=0;k<3;k++)
{
dndt[k]=0;
}
value_type p1,p2,g1,g2,g3,g4,Tp,Tx;
state_type Kt(jmax, 0.0);
int i=0;
for (int j=0;j<jmax;j++)
{

 p1=Tab[0][j];
 p2=Tab[1][j];
 g1=Tab[2][j];
 g2=Tab[3][j];
 g3=Tab[4][j];
 g4=Tab[5][j]; 
 Tp=(p2==10 or g3==10)?Te:Tg;

 Kt[j]={Tab[6][j]*pow(Tp,Tab[7][j])*exp(-Tab[8][j]/Tp)};
/*if (j==0)
{
cerr<<Kt[j]<<endl;

} verification k ok !*/

if (p1==200)
{
Tx=n_Ar*n[p2]*Kt[j];
}
else if (p2==100)
{
Tx=n[p1]*Kt[j];
}
else if (p2==10)
{
Tx=n[p1]*n_e*Kt[j];
}
else
{
Tx=n[p1]*n[p2]*Kt[j];
}


if(p1!=200) {dndt[p1]=dndt[p1]-Tx;}
if(p2!=100 and p2!=10) {dndt[p2]=dndt[p2]-Tx;}
if(g1!=200) {dndt[g1]=dndt[g1]+Tx;}
if(g2!=100) {dndt[g2]=dndt[g2]+Tx;}
if(g3!=100 and g3!=10) {dndt[g3]=dndt[g3]+Tx;}
if(g4!=100) {dndt[g4]=dndt[g4]+Tx;}


}
}
};

/*struct jacobian
{
    void operator()( const state_type &n , matrix_type &jacobi , const value_type &t , state_type &dfdt ) const
    { 
        jacobi( 0 , 0 ) = 0;
        jacobi( 0 , 1 ) = k6*n_e;
        jacobi( 0 , 2 ) = 0;
        jacobi( 1 , 0 ) = 0;
        jacobi( 1 , 1 ) = -k6*n_e;
        jacobi( 1 , 2 ) = 0;
        jacobi( 2 , 0 ) = 0;
        jacobi( 2 , 1 ) = k6*n_e;
        jacobi( 2 , 2 ) = 0;
        dfdt( 0 ) = 0.0;
        dfdt( 1 ) = 0.0;
        dfdt( 2 ) = 0.0;


   }

};*/




void write_density( const value_type t, const value_type Te, const state_type &n)
{
    cout << t << '\t' <<Te <<'\t' 
              << n[0] << '\t' << n[1]<<'\t'<<n[2] << '\t'
              << endl;
}

int main(int argc, char **argv)
{
ifstream fichier_k ("/home/cacot/Documents/rosenbrock/auto.dat");

Tab = new double*[imax];
Tab[0] = new double[imax*jmax];

	for(int i=1;i<imax;i++)
	{
	Tab[i]=Tab[i-1]+jmax;
	}

if(fichier_k)
{
    //Tout est prêt pour la lecture.
cerr<<"fichier ouvert"<<endl;


for(int j=0;j<jmax;j++)
{

           	fichier_k>>Tab[0][j]>>Tab[1][j]>>Tab[2][j]>>Tab[3][j]>>Tab[4][j]
		>>Tab[5][j]>>Tab[6][j]>>Tab[7][j]>>Tab[8][j];
       		
}
fichier_k.close();

}
else
{
    cerr << "ERREUR: Impossible d'ouvrir le fichier en lecture." << endl;
}
cout <<"t"<<'\t'<<"Te"<<'\t'<<"SiH3"<<'\t'<<"SiH4"<<'\t'<<"H"<<endl;

// Time variables
  value_type t = 0.0;
  value_type dt = 1.0e-8;
  value_type Tmax = 20.e-3;
  value_type NT = Tmax/dt;


  // initial values
  value_type Te = 3.0;

  // Density vectors and initial condition
  state_type n_ini(Nbr_espece, 0.0); // initial conditions
  n_ini[0] = 0;
  n_ini[1] = n_SiH4_ini;
  n_ini[2]=0;
 
  state_type n_new(Nbr_espece, 0.0);  // first step same as initial conditions
  n_new = n_ini;
  state_type n_err(Nbr_espece, 0.0); //error

 

  cerr << "\n[ii] Initial Temperature  = " << Te << endl;


  stiff_system* ssys = new stiff_system();

  auto stepper = make_controlled(1.0e-8, 1.0e-8,
                                runge_kutta_cash_karp54< state_type >());
  
//cerr<<"patate 1"<<endl;

size_t num_of_steps = integrate_adaptive(stepper, *ssys, n_new,
                                           t, Tmax, 0.01,
 cout << phoenix::arg_names::arg2 << '\t'
                             << phoenix::arg_names::arg1[0] << '\t'
                             << phoenix::arg_names::arg1[1] << '\t'
		 << phoenix::arg_names::arg1[2] <<'\t'						
                             << phoenix::arg_names::arg1[1]
                            *phoenix::arg_names::arg1[1] << "\n" );
  
//cerr<<"patate 2"<<endl;




 cerr << "\n[ii] Number of steps: " <<  num_of_steps << endl;

//cerr << "SiH3"<<'\t'<<n_new[0]<<'\n'<<"SiH4"<<'\t'<<n_new[1]<<'\n'<<"H"<<'\t'<<n_new[2]<<endl;
//float Si=(n_new[0]+n_new[1])/n_SiH4_ini;

//cerr<<"Si="<<Si<<endl;

//float H=(3*n_new[0]+4*n_new[1]+n_new[2])/(4*n_SiH4_ini);

//cerr<<"H="<<H<<endl;

// deallocate memory for ssys
  delete(ssys);

  return 0;

}  
