#include <stdio.h>
#include <math.h>


#include "glfem.h"
// Théo Denis 27411800
double integrate(double x[3], double y[3], double (*f) (double, double))
{
    double I = 0;
    double xLoc[3];
    double yLoc[3];

    // données initiales pour sommets (0,0) (1,0) (0,1)
    double Xi[3] = {1.0/6.0, 1.0/6.0, 2.0/3.0};
    double Eta[3] = {1.0/6.0, 2.0/3.0, 1.0/6.0};

    // calcul du jacobien dont la valeur absolue est le rapport entre aire quelconque et aire triangle normal
    double J = fabs((x[1]-x[0])*(y[2]-y[0]) - (y[1]-y[0])*(x[2]-x[0]));

    // calcul des nouveaux sommets
    for(int i=0; i<3; i++){
        xLoc[i] = x[0]*(1-Xi[i]-Eta[i]) + x[1]*Xi[i] + x[2]*Eta[i];
        yLoc[i] = y[0]*(1-Xi[i]-Eta[i]) + y[1]*Xi[i] + y[2]*Eta[i];
    }
    
    // calcul de l'intégral
    for(int j=0; j<3; j++){
        I += (1.0/6.0)*f(xLoc[j], yLoc[j]);
    }

    // application du jacobien pour la transformation de domaine
    I = I*J;

    // Pour dessiner l'element, les sommets du triangle :-)
    // Decommenter la ligne pour dessiner aussi les points d'integration
    //

    glfemSetColor(GLFEM_BLACK); glfemDrawElement(x,y,3);
    glfemSetColor(GLFEM_BLUE);  glfemDrawNodes(x,y,3);
    glfemSetColor(GLFEM_RED);   glfemDrawNodes(xLoc,yLoc,3);


    return I;
}

double integrateRecursive(double x[3], double y[3], double (*f)(double,double), int n)
{
    double I = 0;
    if(n>0){
        double x_0[4][3] = {{0.0, 0.0, 0.5}, {0.0, 0.0, 0.5}, {0.5, 0.5, 1.0}, {0.0, 0.5, 0.5}};
        double y_0[4][3] = {{1.0, 0.5, 0.5}, {0.5, 0.0, 0.0}, {0.5, 0.0, 0.0}, {0.5, 0.0, 0.5}};

        double xLoc[4][3];
        double yLoc[4][3];
        
        for(int i=0; i<4; i++){
            for(int j=0; j<3; j++){
                xLoc[i][j] = x[0]*(1-x_0[i][j]-y_0[i][j]) + x[1]*x_0[i][j] + x[2]*y_0[i][j];
                yLoc[i][j] = y[0]*(1-x_0[i][j]-y_0[i][j]) + y[1]*x_0[i][j] + y[2]*y_0[i][j];
            }
            I += integrateRecursive(xLoc[i], yLoc[i], f, n-1);
        }
    } else{
        I = integrate(x,y,f);
    } 
        
    return I;
}