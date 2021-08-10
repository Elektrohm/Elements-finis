# include "fem.h"
// Théo Denis 27411800

femProblem *advdiffNew(double epsilon, double beta, double zeta, int size, int sizePlot)
{
    femProblem *myProblem = malloc(sizeof(femProblem));
    myProblem->epsilon = epsilon;
    myProblem->beta    = beta;
    myProblem->zeta    = zeta;

    myProblem->size = size;    
    myProblem->X = malloc(sizeof(double) * size); 
    myProblem->U = malloc(sizeof(double) * size);
    
    myProblem->sizePlot = sizePlot;    
    myProblem->x = malloc(sizeof(double) * sizePlot); 
    myProblem->u = malloc(sizeof(double) * sizePlot); 
    myProblem->uregime = malloc(sizeof(double) * sizePlot); 
 
    advdiffReset(myProblem);


    return myProblem;
}

void advdiffReset(femProblem *myProblem)
{
    double h;
    int i;

    h = 1.0/(myProblem->size-1);
      for (i=0; i < myProblem->size; i++) {
        myProblem->X[i] = i*h; 
        myProblem->U[i] = 0.0; }
    myProblem->U[0] = 1.0;


    h = 1.0/(myProblem->sizePlot-1);
      for (i=0; i < myProblem->sizePlot; i++) {
        myProblem->x[i] = i*h; 
        myProblem->u[i] = 0.0;        
        myProblem->uregime[i] = 0.0; }
}

void advdiffFree(femProblem *myProblem)
{
    free(myProblem->X);
    free(myProblem->U);
    free(myProblem->x);
    free(myProblem->u);
    free(myProblem->uregime);
    free(myProblem);
}

void advdiffSteadySolution(femProblem *myProblem)
{
    int  size = myProblem->sizePlot, i;
    double *u = myProblem->uregime;
    double *x = myProblem->x;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double pe  = beta/epsilon; 
  
    
    if (beta != 0.0) 
       for (i=0; i < size; i++) 
          u[i] = (exp(pe) - exp(pe*x[i]) )/(exp(pe) - 1.0); 
    else 
       for (i=0; i < size; i++) 
          u[i] = 1 - x[i];          

}

# ifndef NOANALYTIC

void advdiffSolution(femProblem *myProblem, double time)
{
    int  size = myProblem->sizePlot, i,j;
    double *u = myProblem->u;
    double *x = myProblem->x;
    double epsilon = myProblem->epsilon;
    double beta = myProblem->beta;
    double be = beta/epsilon;

    double *uregime = myProblem->uregime;
    advdiffSteadySolution(myProblem);

    // ajoute la solution de régime
    for (i=0; i < size; i++) 
        u[i] = uregime[i];

    // ajoute la solution transitoire
    for (j=1; j < 500; j++) {
        double k = M_PI * j; 
        double a = be/2;
        double b = (beta*beta)/(4*epsilon) + epsilon*k*k;
        double coeff = (8*k)/(4*k*k+be*be);
        for (i=0; i < size; i++) {
                u[i] = u[i] - coeff * sin(k*x[i]) * exp(a*x[i]-b*time); }}
}

# endif
# ifndef NOTIMESTEP

double advdiffComputeTimeStep(femProblem *myProblem)
{
    int size = myProblem->size;
    double beta    = myProblem->beta;
    double epsilon = myProblem->epsilon;
    double zeta = myProblem->zeta ;

    double h = 1.0/(size-1);
    double dt = fmin( (zeta*h*beta + 2*epsilon)/(beta*beta), (h*h)/(zeta*h*beta + 2*epsilon)); 
    
    return 1.35*dt;
}

# endif
# ifndef NORIGHTHANDSIDE

void advdiffComputeRightHandSide(femProblem *myProblem, double *U, double *F)
{
    int  size = myProblem->size, i;
    double epsilon = myProblem->epsilon;
    double beta = myProblem->beta;
    double zeta = myProblem->zeta;
    double h = 1.0/(size-1);
  
    F[0] = 0;
    for (i=1; i < size-1; i++) 
        F[i] = (zeta-1)*beta*(U[i+1]-U[i-1])/(2*h) +  zeta*beta*(U[i-1]-U[i])/h + epsilon*(U[i+1]-2*U[i]+U[i-1])/(h*h);
    F[size-1] = 0;
}

# endif
# ifndef NOEULER

void advdiffUpdateEuler(femProblem *myProblem, double dt)
{
    int size = myProblem->size;
    int i;
    double *U = myProblem->U;
    double *F = malloc(sizeof(double)*size);

    advdiffComputeRightHandSide(myProblem, U, F);

    for(i=0; i<size; i++){
        U[i] += dt*F[i];
    }

    free(F);
}

# endif
# ifndef NORK

void advdiffUpdateRungeKutta(femProblem *myProblem, double dt)
{
 
    int size = myProblem->size;
    int i,k;
    double *U = myProblem->U;
    double *Uold = malloc(sizeof(double)*size);
    double *Uk = malloc(sizeof(double)*size);
    double *K = malloc(sizeof(double)*size);

    // coefficient devant K et pour f(U+coeff*dt*K)
    const double Kcoef[4] = {1.0/6.0, 2.0/6.0, 2.0/6.0, 1.0/6.0};
    const double dtcoef[4] = {0.0, 0.5, 0.5, 1.0};

    for(i=0; i<size; i++){
        // K1 = 0 par sécurité, je sais pas si 0*random = 0 pour 1ère exécu k=0
        K[i] = 0.0;
        Uold[i] = U[i];
    }

    for(k=0; k<4; k++){
        for(i=0; i<size; i++){
            Uk[i] = Uold[i] + dt*dtcoef[k]*K[i];
        }
        // récupère K 
        advdiffComputeRightHandSide(myProblem,Uk,K);
        
        // somme sur U
        for(i=0; i<size; i++){
            U[i] += dt*Kcoef[k]*K[i];
        }
    }

    free(K);
    free(Uk);
    free(Uold);
}

# endif