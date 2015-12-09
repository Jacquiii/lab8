#include <iostream>
#include <cmath>

using namespace std; 

// f(y)
void f(double *f, const double*  const y){
  f[0]=y[1];
  f[1]=-y[0]/sqrt(1+ y[0]*y[0]);
}
   
 

int main(){

const double dt=0.01;
const double T= 20;
const int N= T/dt;

//const double p0=5;




double k1[2];
double k2[2];
double k3[2];
double k4[2];
double ytemp[2];
double yn;


for(double p0 = 0; p0<5; p0+=0.1){
  double y[2];
  y[0] = p0;
  y[1] = 0;
  int counter=0; //Zählt Durchgänge
  
  for(int i=0;i<N;i++){
    
    f(k1,y);  // k1 = f(y_n)
    
    ytemp[0] = y[0] + dt/2 * k1[0];
    ytemp[1] = y[1] + dt/2 * k1[1];

    f(k2,ytemp); // k2 = f( y_n + dt/2 * k1)
      
    ytemp[0] = y[0] + dt/2 * k2[0];
    ytemp[1] = y[1] + dt/2 * k2[1];

    
    f(k3,ytemp); // k3 = f(y_n + dt/2*k2);

    ytemp[0] = y[0] + dt * k3[0];
    ytemp[1] = y[1] + dt * k3[1];
  

    f(k4, ytemp); // k4 = f(y_n + dt * k3)
    yn=y[0];  
    y[0] = y[0] + dt/6 * (k1[0]+2*k2[0]+2*k3[0]+k4[0]);
    y[1] = y[1] + dt/6 * (k1[1]+2*k2[1]+2*k3[1]+k4[1]);
    
    
    if(yn*y[0]<0 && counter ==0){
      //Püfe auf VZW
      double ThetaL, ThetaR, ThetaM, I=1.0;
      ThetaR=1; ThetaL = 0;
      //Bisection
      while(I>1e-10){
	ThetaM = (ThetaR+ThetaR)/2;
	
	double b1 = ThetaM  - 3.0/2*ThetaM*ThetaM + 2.0/3*ThetaM*ThetaM*ThetaM;
	double b2 = ThetaM*ThetaM -2.0/3*ThetaM*ThetaM*ThetaM;
	double b3 = b2;
	double b4 = ThetaM*ThetaM/2 +2.0/3*ThetaM*ThetaM*ThetaM;
	
	I = yn + dt*(k1[0]*b1+k2[0]*b2+k3[0]*b3+k4[0]*b4);
	
	if(I*yn>0)	{ThetaL=ThetaM; yn = I;}
	if(I*y[0]>0)	{ThetaR=ThetaM; y[0]=I;}
	}  
      cout << p0 << "\t" << (i+ThetaM)*4*dt << endl;
      counter++;
    }
   //cout << (i+1)*dt <<"\t"<< y[0]<< "\t" << y[1] <<endl;
  }
}
 return 0;
}
