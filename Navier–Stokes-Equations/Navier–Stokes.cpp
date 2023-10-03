#include <iostream>
#include <stdio.h>
#include <stdlib.h> 
#include <time.h>   
#include <math.h>   
#include <fstream>

int main()
{
    int n, i, j, d;

    n = 200;                                                    // Defining the grid-size

    double u[n][n + 1], un[n][n + 1], uc[n][n + 1];             // Initializing all needed arrays
    double v[n + 1][n], vn[n + 1][n], vc[n + 1][n];             // u for horizontal velocity
    double p[n + 1][n + 1], pn[n + 1][n + 1], pc[n + 1][n + 1]; // v for vertical velocity
    double m[n + 1][n + 1];                                     // p for vorticity
    double dx, dy, dt, delta, err, Re, t;
    int dmax = 100000;

    std::ofstream myfile("Navier_stokes.dat");                  // Opening data file for later use.
    std::ofstream myfile2("Navier_stokes2.dat");                // Opening data file for later use.
    std::ofstream myfile3("Navier_stokes3.dat");
    std::ofstream myfile4("Navier_stokes4.dat");                // Opening data file for later use.
    myfile.precision(17);                                       // Defining precision for each file.
    myfile2.precision(17);
    myfile3.precision(17);
    myfile4.precision(17);
    std::string buf;                                            // String stream buffer.

    d = 1;                                                      // intializing step counter, d
    dx = 1.0 / (n - 1);                                         // dx = 1/gridsize
    dy = dx;                                                    // we are using a square, therefore dy=dx
    dt = 0.00002;                                               // dt is chosen to fit the gridsize and a comfortable divergence number/s
    delta = 4.0;                                                // WHAT IS THIS I DO NOT KNOW
    err = 2.0;                                                  // Setting intial error above tolerance to start the loop
    Re = 100000.0;                                              // Reynolds number for water: ~32e3
    double Divergence1, Divergence2;                            //
    Divergence1 = dt / dx;
    Divergence2 = (1.0 / Re) * dt / (dx * dx);
    printf("Divergence numbers are %lf and %lf\n", Divergence1, Divergence2);

    // ---------------------------------------------------------------------------------------
    // u,v,p startup and setting Initial Conditions
    //
    for(i = 0; i <= (n - 1); i++)
    {
        for(j = 0; j <= (n); j++)
        {

            u[i][j] = 0.0;
            u[i][n] = 1.0;
            u[i][n - 1] = 1.0;
        }
    }

    for(i = 0; i <= (n); i++)
    {
        for (j = 0; j <= (n - 1); j++)
        {
            v[i][j] = 0.0;
        }
    }

    for (i = 0; i <= (n); i++)
    {

        for (j = 0; j <= (n); j++)
        {
            p[i][j] = 1.0;
        }
    }

    t = 0.0;
    while(d < 8000000)
    {
        // Interior Points Calculation
        for(i = 1; i <= (n - 2); i++)
        {
            for (j = 1; j <= (n - 1); j++)
            {
                un[i][j] = u[i][j] - dt * ((u[i + 1][j] * u[i + 1][j] - u[i - 1][j] * u[i - 1][j]) / 2.0 / dx + 0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i][j] + u[i][j - 1]) * (v[i + 1][j - 1] + v[i][j - 1])) / dy) - dt / dx * (p[i + 1][j] - p[i][j]) + dt * 1.0 / Re * ((u[i + 1][j] - 2.0 * u[i][j] + u[i - 1][j]) / dx / dx + (u[i][j + 1] - 2.0 * u[i][j] + u[i][j - 1]) / dy / dy);
            }
        }

        // B.C.
        for(j = 1; j <= (n - 1); j++)
        {
            un[0][j] = 0.0;
            un[n - 1][j] = 0.0;
        }
        // B.C.
        for(i = 0; i <= (n - 1); i++)
        {
            un[i][0] = -un[i][1];
            un[i][n] = 2.0 - un[i][n - 1];
        }

        // Solves v-momentum
        for(i = 1; i <= (n - 1); i++)
        {
            for(j = 1; j <= (n - 2); j++)
            {
                vn[i][j] = v[i][j] - dt * (0.25 * ((u[i][j] + u[i][j + 1]) * (v[i][j] + v[i + 1][j]) - (u[i - 1][j] + u[i - 1][j + 1]) * (v[i][j] + v[i - 1][j])) / dx + (v[i][j + 1] * v[i][j + 1] - v[i][j - 1] * v[i][j - 1]) / 2.0 / dy) - dt / dy * (p[i][j + 1] - p[i][j]) + dt * 1.0 / Re * ((v[i + 1][j] - 2.0 * v[i][j] + v[i - 1][j]) / dx / dx + (v[i][j + 1] - 2.0 * v[i][j] + v[i][j - 1]) / dy / dy);
            }
        }

        // B.C.
        for(j = 1; j <= (n - 2); j++)
        {
            vn[0][j] = -vn[1][j];
            vn[n][j] = -vn[n - 1][j];
        }

        // B.C.
        for(i = 0; i <= (n); i++)
        {
            vn[i][0] = 0.0;
            vn[i][n - 1] = 0.0;
        }

        // ---------------------------------------------------------------------------------------
        // continuity equation
        for(i = 1; i <= (n - 1); i++)
        {
            for(j = 1; j <= (n - 1); j++)
            {
                pn[i][j] = p[i][j] - dt * delta * ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
            }
        }

        // B.C.
        for(i = 0; i <= (n); i++)
        {
            pn[i][0] = pn[i][1];
            pn[i][n] = pn[i][n - 1];
        }

        for(j = 0; j <= (n); j++)
        {
            pn[0][j] = pn[1][j];
            un[n][j] = pn[n - 1][j];
        }

        // ---------------------------------------------------------------------------------------
        err = 0.0;
        for(i = 1; i <= (n - 1); i++)
        {
            for(j = 1; j <= (n - 1); j++)
            {
                m[i][j] = ((un[i][j] - un[i - 1][j]) / dx + (vn[i][j] - vn[i][j - 1]) / dy);
                err = err + fabs(m[i][j]);
            }
        }

        // residual[step] = log10(error);
        if(d % 2000 == 1)
        {
            printf("Error is %5.5lf, time is %5.5lf of 160 seconds for the step %d\n", err, t, d);
        }

        for(i = 0; i <= (n - 1); i++)
        {
            for(j = 0; j <= (n); j++)
            {
                u[i][j] = un[i][j];
            }
        }

        for(i = 0; i <= (n); i++)
        {
            for(j = 0; j <= (n - 1); j++)
            {
                v[i][j] = vn[i][j];
            }
        }

        for(i = 0; i <= (n); i++)
        {
            for(j = 0; j <= (n); j++)
            {
                p[i][j] = pn[i][j];
            }
        }

        d = d + 1;
        t = t + dt;


        if(d % 2000 == 1)
        {
            for(i = 1; i <= (n - 1); i++)
            {
                for(j = 1; j <= (n - 1); j++)
                {
                    myfile << p[i][j] << '\n';
                }
            }

            for(i = 1; i <= (n - 1); i++)
            {
                for(j = 1; j <= (n - 1); j++)

                    myfile2 << u[i][j] << '\n';
            }

            for(i = 1; i <= (n - 1); i++)
            {
                for(j = 1; j <= (n - 1); j++)

                    myfile3 << v[i][j] << '\n';
            }

            myfile4 << t << '\n';
        }
    }

    myfile.close();
}
