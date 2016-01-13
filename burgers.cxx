#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <cmath>
//---------------------------------------
using namespace std;
//---------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N);
void initialize(double* const u1, double* const u0, const double dx,const double dt, const double xmin,
                const int N);
void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N);

//---------------------------------------
int main(){

  const double tEnd = 0.4 ;


  const int N  = 64;
  const double xmin = 0;
  const double xmax = 1;
  const double dx = (xmax-xmin)/N ;
  double dt = 0.001;
  double t = 0;
  const int Na = 10;
  const int Nk = int(tEnd/Na/dt);

  double* u0 = new double[N];
  double* u1 = new double[N];
  double* u2 = new double[N];
  double* h;
  stringstream strm;

  initialize(u1,u0,dx,dt, xmin,N);

  writeToFile(u0, "u_0", dx, xmin, N);

  cout << "Nk = " << Nk << endl;
  double komisch,u;
   double x ;
  
  ofstream result("result.dat"); // result = analytisch
 for(int i=0; i < N; i++)
 {
   x = xmin + i*dx;
   u = sin(2.0*M_PI*x);
   komisch = x + u*tEnd;
   result << komisch << '\t' << u << endl;
 }
  result.close();

  for(int i=1; i<=Na; i++)
  {
   for(int j=0; j<Nk; j++){

      step(u2, u1, u0, dt, dx, N);     // step + swap here
      h = u0;
      u0 = u1;
      u1 = u2;
      u2=h;
      
      
      t +=dt;
   }
   strm.str("");
   strm << "u_" << i;
   writeToFile(u0, strm.str(), dx, xmin, N);
  }

  cout << "t = " << t << endl;

  delete[] u0;
  delete[] u1;
  delete[] u2;
  return 0;
}
//-----------------------------------------------
void step(double* const u2, const double* const u1,const double* const u0,
          const double dt, const double dx, const int N)
{

 u2[0]= u0[0]-dt*u1[0]/dx*(u1[+1]-u1[N-1]) ; 
  
  for(int i=1 ; i<N-1 ; i++) {
  u2[i]=u0[i]-dt*u1[i]/dx*(u1[i+1]-u1[i-1]) ; 

  }
  
  u2[N-1]=u0[N-1]-dt*u1[N-1]/dx*(u1[0]-u1[N-2]) ; // der N-te Wert ist der 0.te Wert, da wir von 0-63 gehen und N=64 (0.te)
}
//-----------------------------------------------
void initialize(double* const u1, double* const u0, const double dx,
                const double dt, const double xmin,  const int N)
{
  
  //Anfangsbed.
   double u,ux, uxx;
   for(int i=0; i<N; i++)
   {
     double x = xmin + i*dx;
     
     u = sin(2.0* M_PI*x);
     ux = 2.0*M_PI*cos(2.0* M_PI*x);
     uxx = -4.0*M_PI*M_PI*sin(2.0* M_PI*x);

     u1[i] = u; 
     u0[i] = u1[i] + dt*u1[i]*ux+0.5*dt*dt*(u1[i]*(2.0*ux*ux)+u1[i]*uxx);
   }
}
//-----------------------------------------------
void writeToFile(const double* const u, const string s, const double dx,
                 const double xmin, const int N)
{
   ofstream out(s.c_str());
   for(int i=0; i<N; i++){
     double x = xmin + i * dx;
     out << x << "\t" << u[i] << endl;
   }
   out.close();
}
