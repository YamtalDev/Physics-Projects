//----------------------------------3D Diffusion equation c++--------------------------------------------------------------------
// This code solves the three-dimensional diffusion equation using numerical tools for the patron of differential equations and numerical derivatives.
// The domain where the diffusiom occur is in 1x1x1 cubical box, Begining by defining the time and place interval, and the stability corient.
// Set initial conditions for concentration through all the domain (i,j,k) at time t = 0 to be -0- to initialize diffusion.
// Creating temporery arrey and initializing it to. Defining pollution at some point in the domain.
// Creating the numeric derivatives to itterat over them using nested loops.
// Printing the results to a data file for later use.
// Finishing with the loop over time.
//-------------------------------------------------------------------------------------------------------------------------------

#include <iostream>                  // Including Standard Input / Output Streams Library.
#include <fstream>                   // Including Input/output file stream class iostream.
#include <cstdlib>                   // Including several general purpose functions stdlib.
#include <string>                    // Including string to represent sequences of characters.
#include <cmath>                     // Including cmath functions to compute common mathematical operations.

using namespace std;

int main(){                          // Opening int main to preform calculations.

//-------------------------------------- Decleration of variables ----------------------------------------------------------------
long double pi = 3.14592654;
const int P = 55;                    // Size of wanted arrey.
long double C[P][P][P];              // Creating the inital arrey of the diffusion.
long double C_t[P][P][P];            // Creating a temporery arrey to itterat over time on it.
int Ny,Nx,Nz,Nt;                     // Declering the the number of steps in every axis and time.
int i,j,k,l;                         // Declering the loop integers.
long double dx,dy,dz,dt;             // Differentials of x,y,z,t axis.
long double Sx,Sy,Sz;
long double Cx,Cz,Cy;
long double x,y,z,t;                 // 3D coardinates.
long double u,v,w;                   // Currents speed in the u,v,w coardinates.
long double coreant,coreant2, ana;   // Stability coreants.
long double Difx,Dify,Difz;          // Diffusibility constants.

//-------------------------------------- Variables calculation --------------------------------------------------------------------

ofstream  myfile("diff.dat");        // Opening data file for later use.

myfile.precision(17);                // Defining our precision.

string buf;                          // String stream buffer.

cout <<" Enter final time:"<<endl;   // Asking the user for the final time.

cin >>t;getline(cin,buf);            // Entering the value of time to the memory.

Nx=20; Ny=20; Nz=20; Nt=2001;        // Creating number of steps in the 3D axis and precision.

Difx= 0.002;                       // The Diffusibility constant in the x direction.

Dify= 0.002;                       // The Diffusibility constant in the y direction.

Difz= 0.002;                        // The Diffusibility constant in the z direction.

dx=1.0/(Nx-1);                       // Calculating axis Differentials.

dy=1.0/(Ny-1);                       // Calculating axis Differentials.

dz=1.0/(Nz-1);                       // Calculating axis Differentials.

dt=t/(Nt-1);                         // Calculating time Differentials.

u=0.0; v=0.0; w=0.0;                 // Currents speeds in the u,v,w grid.

Sx=(Difx*dt)/(pow(dx,2));            // Creating the constant parameters of the numeric calculation.

Sy=(Dify*dt)/(pow(dy,2));            // Creating the constant parameters of the numeric calculation.

Sz=(Difz*dt)/(pow(dz,2));            // Creating the constant parameters of the numeric calculation.

Cx=(u*dt)/(dx);                      // Creating the constant parameters of the numeric calculation.

Cy=(v*dt)/(dy);                      // Creating the constant parameters of the numeric calculation.

Cz=(w*dt)/(dz);                      // Creating the constant parameters of the numeric calculation.


coreant=Sx+Sy+Sz;                                                            // Calculating coreant for stability.

coreant2 = (pow(Cx,2)/Sx)+(pow(Cy,2)/Sy)+(pow(Cz,2)/Sz);                     // Calculating second coreant for stability.

if(coreant2 > 3.0){cout << "cnt2 =" << coreant2 << endl; return(0);}         // If statments to make sure coreant will be less then 3 for good stability.

if(coreant > 0.5){cout << "Cnt=" << coreant << endl;return(0);}              // If statments to make sure coreant will be less then 0.5 for good stability.

cout <<"  dt=  "<<dt<<"  dx= "<<dx<<"  dy= "<<dy<<"  dz= "<<dz<< endl;       // Showing the user the Differentials.

cout << "cnt2 = " <<coreant2<<" cnt = " <<coreant<< endl;                    // Showing the user the 2 coreants.

//-------------------------------------- Inital conditions and bounderies ---------------------------------------------------------//

                                     // Opening the nested for loops to calculat the 3D grid.

for(k=1;k<=Nz-1;k++)                 // Loop for z axis.

   for(i=1;i<=Nx-1;i++)              // Loop for x axis.

      for(j=1;j<=Ny-1;j++){          // Loop for y axis.



C[i][j][k]=0.0;                      // Inital condition for C(0,i,j,k)=0.0 to preset all the grid.

C_t[i][j][k]=0.0;                    // Doing the same with the temporery arrey.


}



l=0;
C[10][10][10]=300;
do{


    for(i=1;i<Nx-1;i++)
        for(j=1;j<Nx-1;j++)
            for(k=1;k<Nz-1;k++){










C_t[i][j][k]=(Sx+(Cx/2.0))*C[i-1][j][k]+(Sy+(Cy/2.0))*C[i][j-1][k]
+(Sz+(Cz/2.0))*C[i][j][k-1]+(Sx-(Cx/2.0))*C[i+1][j][k]
+(Sy-(Cy/2.0))*C[i][j+1][k]+(Sz-(Cz/2.0))*C[i][j][k+1]
+(1.0-2.0*(coreant))*C[i][j][k];

C[i][j][k]=C_t[i][j][k];

//----------------------------------------------------------------- ft bx

C_t[Nx-1][j][k]=C[Nx-1][j][k]
-(Cx/2)*(3.0*C[Nx-1][j][k]-4.0*C[Nx-2][j][k]+3.0*C[Nx-3][i][k])
-(Cy/2)*(C[Nx-1][j+1][k]-C[Nx-1][j-1][k])
-(Cz/2)*(C[Nx-1][j][k+1]-C[Nx-1][j][k-1])
+Sx*(C[Nx-1][j][k]-2.0*C[Nx-2][j][k]+C[Nx-3][j][k])
+Sy*(C[Nx-1][j+1][k]-2.0*C[Nx-1][j][k]+C[Nx-1][j-1][k])
+Sz*(C[Nx-1][j][k+1]-2.0*C[Nx-1][j][k]+C[Nx-1][j][k-1]);

//----------------------------------------------------------------- ft by

C_t[i][Ny-1][k]=C[i][Ny-1][k]
-(Cx/2)*(C[i+1][Ny-1][k]-C[i-1][Ny-1][k])
-(Cy/2)*(3.0*C[i][Ny-1][k]-4.0*C[i][Ny-2][k]+3.0*C[i][Ny-2][k])
-(Cz/2)*(C[i][Ny-1][k+1]-C[i][Ny-1][k-1])
+Sx*(C[i+1][Ny-1][k]-2.0*C[i][Ny-1][k]+C[i-1][Ny-1][k])
+Sy*(C[i][Ny-1][k]-2.0*C[i][Ny-2][k]+C[i][Ny-3][k])
+Sz*(C[i][Ny-1][k+1]-2.0*C[i][Ny-1][k]+C[i][Ny-1][k-1]);

//---------------------------------------------------------------- ft bz


C_t[i][j][Nz-1]=C[i][j][Nz-1]
-(Cx/2)*(C[i+1][j][Nz-1]-C[i-1][j][Nz-1])
-(Cy/2)*(C[i][j+1][Nz-1]-C[i][j-1][Nz-1])
-(Cz/2)*(3.0*C[i][j][Nz-1]-4.0*C[i][j][Nz-2]+3.0*C[i][j][Nz-3])
+Sx*(C[i+1][j][Nz-1]-2.0*C[i][j][Nz-1]+C[i-1][j][Nz-1])
+Sy*(C[i][j+1][Nz-1]-2.0*C[i][j][Nz-1]+C[i][j-1][Nz-1])
+Sz*(C[i][j][Nz-1]-2.0*C[i][j][Nz-2]+C[i][j][Nz-3]);

//----------------------------------------------------------------- ft fx

C_t[1][j][k]=C[1][j][k]
-(Cx/2)*(C[1][j][k]-4.0*C[2][j][k]+3.0*C[3][i][k])
-(Cy/2)*(C[1][j+1][k]-C[1][j-1][k])
-(Cz/2)*(C[1][j][k+1]-C[1][j][k-1])
+Sx*(C[1][j][k]-2.0*C[2][j][k]+C[3][j][k])
+Sy*(C[1][j+1][k]-2.0*C[1][j][k]+C[1][j-1][k])
+Sz*(C[1][j][k+1]-2.0*C[1][j][k]+C[1][j][k-1]);

//----------------------------------------------------------------- ft fy

C_t[i][1][k]=C[i][1][k]
-(Cx/2)*(C[i+1][1][k]-C[i-1][1][k])
-(Cy/2)*(C[i][1][k]-4.0*C[i][2][k]+3.0*C[i][3][k])
-(Cz/2)*(C[i][1][k+1]-C[i][1][k-1])
+Sx*(C[i+1][1][k]-2.0*C[i][1][k]+C[i-1][1][k])
+Sy*(C[i][1][k]-2.0*C[i][2][k]+C[i][3][k])
+Sz*(C[i][1][k+1]-2.0*C[i][1][k]+C[i][1][k-1]);

//---------------------------------------------------------------- ft fz

C_t[i][j][1]=C[i][j][1]
-(Cx/2)*(C[i+1][j][1]-C[i-1][j][1])
-(Cy/2)*(C[i][j+1][1]-C[i][j-1][1])
-(Cz/2)*(C[i][j][1]-4.0*C[i][j][2]+3.0*C[i][j][3])
+Sx*(C[i+1][j][1]-2.0*C[i][j][1]+C[i-1][j][1])
+Sy*(C[i][j+1][1]-2.0*C[i][j][1]+C[i][j-1][1])
+Sz*(C[i][j][1]-2.0*C[i][j][2]+C[i][j][3]);

C[0][j][k]=0.0;
C[i][0][k]=0.0;
C[i][j][0]=0.0;
C[Nx][j][k]=0.0;
C[i][Ny][k]=0.0;
C[i][j][Nz]=0.0;
 }

l++;


for(k=1;k<=Nz-1;k++)
    for(i=1;i<=Nx-1;i++)
        for(j=1;j<=Nx-1;j++){


myfile  << C[i][j][k] << '\n';
        }
}while(l<Nt-1);



myfile.close();
}
