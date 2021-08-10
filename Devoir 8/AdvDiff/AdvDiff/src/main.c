/*
 *  main.c
 *  Library for LEPL1110 : Finite Elements for dummies
 *
 *  Copyright (C) 2021 UCL-EPL : Vincent Legat
 *  All rights reserved.
 *
 */

# include "glfem.h"

int main(void)
{   

// 
//
// Diffusion pure      :  epsilon = 0.01  
// Advection-diffusion :  epsilon = 0.01  beta = 1
//     Difference centre -> zeta = 0.0
//     Pure upwind       -> zeta = 1.0
//

    int sizePlot   = 101;   // Discretisation pour les solutions analytiques (le plot :-) 
    int size       = 21;    // Discretisation différences finies
    double beta    = 0.2;
    double epsilon = 0.01;
    
    double dx = 1.0/(size-1);
    double peclet = beta*dx/epsilon;
    double zeta = 1.0/tanh(peclet/2) - 2.0/peclet;
    femProblem *theProblem = advdiffNew(epsilon,beta,zeta,size,sizePlot);
     
    double 	theDiscreteTime = 0.0;
    double 	theStartingTime = 0.0;
    double  theTimeStep;
    int 	  theIteration = 0;
    femTemporalSchemeType theSolver = FEM_EULER;
    theTimeStep = dx*dx/(epsilon*2.00);

//
// Integration hors de l'environnement graphique
// Oui, c'est possible :-)
//
    
	  int i;
    for (i=0; i < 100; i++) {
        advdiffUpdateEuler(theProblem,theTimeStep);
        if (i%10 == 0) printf("Iteration Euler       %2d : %14.7e \n",i,theProblem->U[size/2]); }
    advdiffReset(theProblem);
     
    for (i=0; i < 100; i++) {
        advdiffUpdateRungeKutta(theProblem,theTimeStep);
        if (i%10 == 0) printf("Iteration Runge-Kutta %2d : %14.7e \n",i,theProblem->U[size/2]); }
    advdiffReset(theProblem);
       
//
// Definition de l'environnement graphique
//

    const char theHelpMessage[] = {
    "   [esc] : Exit\n"
    "    R    : Restart and reset zoom, translations \n"
    "    E-K  : Euler - Runge/Kutta scheme \n"
    "    H    : Display or hide keyboard shortcuts \n"};
    printf("\n%s\n",theHelpMessage);
    glfemWindowCreate("EPL1110 : Advection-diffusion equation",480,480,2,
                                       (double[2]){0.0,1.0},(double[2]){0.0,1.0});
    glfemWindowSetHelpMessage(theHelpMessage);                               
                                       
    do 
    {     
        glfemWindowUpdate();
        char theMessage[256];
        char action = glfemGetAction(); 
        
//
//  Gestion de l'animation en temps reel
//    - "theTime" est le temps courant de l'application (sauf si il a été remis a zero :-)
//    - Le facteur entre le temps réel et le temps de l'application permet de controler la vitesse d'execution
//    - Si necessaire, une nouvelle iteration discrete est calculee...
//      C'est donc ici que se trouve "virtuellement" la boucle sur toutes les iterations temporelles
//
//  Pour figer/ne pas figer le resultat a un temps, decommenter les deux lignes ci-dessous
//
        double theTime = (glfwGetTime() - theStartingTime) * 5;   
        double theStop = 2;
        if (theTime >= theStop) theTime = theStop;
        
        if (action != '0') {
            printf(" ==== Command dectected : %c \n",action);
            if (action == 'K')  theSolver = FEM_RK;
            if (action == 'E')  theSolver = FEM_EULER;
            theStartingTime = glfwGetTime(); 
            theDiscreteTime = 0;
    		    theIteration = 0;
    		    advdiffReset(theProblem); }     
  
        advdiffSteadySolution(theProblem);
        advdiffSolution(theProblem,theTime);
        const char *theSolverName = femTemporalSchemeName(theSolver);
        if (theTime >= theDiscreteTime) {
           theIteration += 1;
           theDiscreteTime += theTimeStep; 
           printf("Iteration %s %2d - %.2f : %14.7e \n",theSolverName,theIteration,theDiscreteTime,theProblem->U[size/2]); 
           if (theSolver == FEM_EULER)   advdiffUpdateEuler(theProblem,theTimeStep); 
           if (theSolver == FEM_RK)      advdiffUpdateRungeKutta(theProblem,theTimeStep); }
        sprintf(theMessage,"%s : time = %.2f iteration = %d",theSolverName,theDiscreteTime,theIteration);

//
// Definition du dessin a realiser au temps "theTime"
//
//   - en rouge : solution analytique courante de l'equation de la chaleur (pas d'advection)
//   - en vert : solution analytique de regime du probleme d'advection-diffusion
//   - en bleu : derniere solution discrete obtenue
//
 
        glfemSetColor(GLFEM_RED);   glfemDrawCurve(theProblem->x,theProblem->u,theProblem->sizePlot);
        glfemSetColor(GLFEM_GREEN); glfemDrawCurve(theProblem->x,theProblem->uregime,theProblem->sizePlot);
        glfemSetColor(GLFEM_BLUE);  glfemDrawCurveDiscrete(theProblem->X,theProblem->U,theProblem->size);
        glfemDrawNodes(theProblem->X,theProblem->U,theProblem->size);       
     	  glfemDrawMessage(theMessage,(double[2]){16.0, 30.0});	
              
    } while(!glfemWindowShouldClose());
        
    glfemWindowFree(); 
    advdiffFree(theProblem);
    exit(EXIT_SUCCESS); 

}

