#include"fem.h"
// Théo Denis 27411800

double convergenceSource(double x, double y)
{

//
// A obtenir avec sympy :-)
//
    double f = 2*(x*(x - 1)*(-400*y*(y - 1)*(5*M_SQRT2*(x + y) - 8)
     + 10*M_SQRT2*(2*y - 1)*(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1) 
     + pow(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1, 2)*atan(10*M_SQRT2*(x + y) - 16))
     + y*(y - 1)*(-400*x*(x - 1)*(5*M_SQRT2*(x + y) - 8) + 10*M_SQRT2*(2*x - 1)*
     (4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1) + pow(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1, 2)
     *atan(10*M_SQRT2*(x + y) - 16)))/pow(4*pow(5*M_SQRT2*(x + y) - 8, 2) + 1, 2);
    return -f; 
}

double convergenceSoluce(double x, double y, double *u)
{
    double a = sqrt(2.0);
    u[0] = x*y*(1-x)*(1-y)*atan(20*(x+y)/a - 16);
    u[1] = y*(y - 1)*(10*a*x*(x - 1)
      + (2*x - 1)*(4*pow(5*a*(x + y) - 8, 2) + 1)*atan(10*a*(x + y) - 16))
      /(4*pow(5*a*(x + y) - 8, 2) + 1);
    u[2] = x*(x - 1)*(10*a*y*(y - 1) 
      + (2*y - 1)*(4*pow(5*a*(x + y) - 8, 2) + 1)*atan(10*a*(x + y) - 16))
      /(4*pow(5*a*(x + y) - 8, 2) + 1);
    return u[0];
}

double convergenceEstimateRate(double *errors, int n, double ratio)
{

     double rate = 0;

     for(int i=0; i<n-1; i++){
       rate += ((log10(errors[i]/errors[i+1]))*(1.0/log10(ratio)));
     }
     rate = rate*(1.0/(n-1));
     return rate;
}


void femDiffusionComputeError(femDiffusionProblem *theProblem, 
                                    double(*soluce)(double,double,double*))
{


//  
  femMesh *theMesh = theProblem->mesh;
  femIntegration *theRule = theProblem->rule;
  femDiscrete *theSpace = theProblem->space;
  
  int n = theProblem->space->n;

  int i,j;
  int map[4], dirclette[4];
  double xloc[4],yloc[4],phi[4],dphidx[4], dphidy[4], dphidxsi[4],dphideta[4], Uloc[4], Utilde[4], poubelle[3];

  theProblem->errorSoluceL2 = 0.0;
  theProblem->errorSoluceH1 = 0.0;
  theProblem->errorInterpolationL2 = 0.0;
  theProblem->errorInterpolationH1 = 0.0;

  // boucle sur les éléments du maillages
  for(i=0; i<theMesh->nElem; i++){
    // récupère valeur de l'élément local
    femDiffusionMeshLocal(theProblem,i,map, dirclette,xloc,yloc,Uloc);

    // récupère valeur de Utilde sur chaque noeud
    for(int t=0 ; t<n; t++){
      soluce(xloc[t],yloc[t],poubelle);
      Utilde[t] = poubelle[0];
    }

    // changement de variable 
      for(j=0; j<theRule->n; j++){
          // récupère paramètre intégration élément parents
          double xsi = theProblem->rule->xsi[j];
          double eta = theProblem->rule->eta[j];
          double weight = theProblem->rule->weight[j];

          // récupère fonction de forme
          femDiscretePhi2(theProblem->space, xsi, eta, phi);
          femDiscreteDphi2(theProblem->space, xsi, eta, dphidxsi, dphideta);

          // aide pour le changement de variable
          double dxdxsi = 0.0;
          double dxdeta = 0.0;
          double dydxsi = 0.0;
          double dydeta = 0.0;

          double xmoy = 0.0;
          double ymoy= 0.0;

          double uh = 0.0;
          double utilde = 0.0;

          // variable changement de var 
          double jac = 0.0;

          for(int k=0; k<n; k++){
              xmoy += xloc[k]*phi[k];
              ymoy += yloc[k]*phi[k];
              uh += Uloc[k]*phi[k];
              utilde += Utilde[k]*phi[k];

              dxdxsi += xloc[k]*dphidxsi[k];
              dxdeta += xloc[k]*dphideta[k];

              dydxsi += yloc[k]*dphidxsi[k];
              dydeta += yloc[k]*dphideta[k];
          }

          jac = fabs(dxdxsi*dydeta-dydxsi*dxdeta);
        
          double duhdx = 0.0;
          double duhdy = 0.0;
          double dutildedx = 0.0;
          double dutildedy = 0.0;
          

          // applique le changement de variable
          for(int l=0; l<n; l++){
              dphidx[l] = (1.0/jac)*(dydeta*dphidxsi[l]- dydxsi*dphideta[l]);
              dphidy[l] = (1.0/jac)*(dxdxsi*dphideta[l]- dxdeta*dphidxsi[l]);

              duhdx += Uloc[l]*dphidx[l];
              duhdy += Uloc[l]*dphidy[l];

              dutildedx += Utilde[l]*dphidx[l];
              dutildedy += Utilde[l]*dphidy[l];
          }

          // calcul de u
          soluce(xmoy,ymoy, poubelle);
          double u = poubelle[0];
          double dudx = poubelle[1];
          double dudy = poubelle[2];
         
          // erreur sur uh
          double erroruh = (u-uh)*(u-uh);
          double errorduhdx = (dudx-duhdx)*(dudx-duhdx);
          double errorduhdy = (dudy-duhdy)*(dudy-duhdy);

          // erreur sur utilde
          double errorutilde = (u-utilde)*(u-utilde);
          double errordutildedx = (dudx-dutildedx)*(dudx-dutildedx);
          double errordutildedy = (dudy-dutildedy)*(dudy-dutildedy);

          theProblem->errorSoluceL2 += erroruh*jac*weight;
          theProblem->errorSoluceH1 += (erroruh + errorduhdx + errorduhdy)*jac*weight;
          theProblem->errorInterpolationL2 += errorutilde * jac * weight;
          theProblem->errorInterpolationH1 += (errorutilde + errordutildedx + errordutildedy)*jac*weight;
      }
  }
  theProblem->errorSoluceL2 = sqrt(theProblem->errorSoluceL2);
  theProblem->errorSoluceH1 = sqrt(theProblem->errorSoluceH1);
  theProblem->errorInterpolationL2 = sqrt(theProblem->errorInterpolationL2);
  theProblem->errorInterpolationH1 = sqrt(theProblem->errorInterpolationH1);
}




femMesh *femMeshCreateBasicSquare(int n) 
{
    femMesh *theMesh = malloc(sizeof(femMesh));

    // Toujours le maillage ?l?mentaire :-)
    // A modifier [begin]
        
    // allocation des vecteurs X,Y
    theMesh->nNode = (n+1)*(n+1);
    theMesh->X = malloc(sizeof(double)*theMesh->nNode);
    theMesh->Y = malloc(sizeof(double)*theMesh->nNode);
    double *X = theMesh->X;
    double *Y = theMesh->Y;
    int i,j;
    double val;
    
    // initialisation des noeuds
    for(i=0; i<n+1; i++){
      val = i*1.0/n;
      for(j=0; j<n+1; j++){
        X[(n+1)*i+j] = val;
        Y[i+j*(n+1)] = val;
      }
    }

    // allocation du tableau d'élément
    theMesh->nElem = 2*n*n;
    theMesh->nLocalNode = 3;
    theMesh->elem = malloc(sizeof(int)*3*theMesh->nElem);
    int *elem = theMesh->elem;
    
    int incr = 0;
    int index = 0;

    // initialisation des éléments
    for(i=0; i<n*n; i++){
      for(j=0; j<2; j++){
        elem[3*index] = i+j*(n+1)+incr;
        elem[3*index+1] = i+n+2+incr;
        elem[3*index+2] = elem[3*index+1]-n-(j+1);
        index++;
      }
      if ((i+1)%n==0) incr++;
    }
    
    // A modifier [end]

    theMesh->number = malloc(sizeof(int)*theMesh->nNode); 
    for (int i = 0; i < theMesh->nNode; i++) 
          theMesh->number[i] = i;     
    return theMesh;
}
