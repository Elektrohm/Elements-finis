#include "glfem.h"
#include <time.h>

int main(void)
{  

    printf("\n\n    V : Plot results \n");
    printf("    S : Spy matrix \n");
    printf("    F-B-I : Full solver - Band solver - Iterative solver \n");
    printf("    X-Y-N : Renumbering along x - along y - No renumbering \n");
    
    femSolverType solverType = FEM_BAND;
    femRenumType  renumType  = FEM_YNUM;
    
    char meshFileName[] = "../data/triangles101.txt";  
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles101.txt") 
    // par :
    // ("..\\data\\triangles101.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !

    
    femDiffusionProblem* theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
    clock_t tic = clock();
    femDiffusionCompute(theProblem);  
    femSolverPrintInfos(theProblem->solver);
    printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
    printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
    fflush(stdout);
    
    glfemWindowCreate("EPL1110 : Band Solver",480,480,theProblem->mesh->nNode,theProblem->mesh->X,theProblem->mesh->Y);


    int option = 1;    
    femSolverType newSolverType = solverType;
    femRenumType  newRenumType  = renumType;
    
    do 
    {
        glfemWindowUpdate();    
        
        int testConvergence,w,h;
        char theMessage[256];
        sprintf(theMessage, "Max : %.4f ",femMax(theProblem->soluce,theProblem->size));
        
          
        if (option == 1) {
            glfemSetColor(GLFEM_BLACK);
            glfemPlotSolution(theProblem->mesh,theProblem->soluce);}
        else {
            glfemPlotSolver(theProblem->solver,theProblem->size);}
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0});  
        
        if (solverType != newSolverType || renumType != newRenumType) { 
            solverType = newSolverType;
            renumType = newRenumType;
            femDiffusionFree(theProblem);
            theProblem = femDiffusionCreate(meshFileName,solverType,renumType);
            clock_t tic = clock();
            do {
                femDiffusionCompute(theProblem);  
                femSolverPrintInfos(theProblem->solver); 
                testConvergence = femSolverConverged(theProblem->solver); }
            while ( testConvergence == 0);
            if (testConvergence == -1)  printf("    Iterative solver stopped afer a maximum number of iterations\n");
            printf("    CPU time : %.2f [sec] \n", (clock() - tic) * 1.0 /CLOCKS_PER_SEC);
            switch (renumType) {
                case FEM_XNUM : printf("    Renumbering along the x-direction\n"); break;
                case FEM_YNUM : printf("    Renumbering along the y-direction\n"); break;
                default : break; }
            printf("    Maximum value : %.4f\n", femMax(theProblem->soluce,theProblem->size));
            fflush(stdout); }           
          
          if (glfemGetKey('V'))  {option = 1; glfemWindowResetSize();}
          if (glfemGetKey('S'))   option = 0;
          if (glfemGetKey('F'))   newSolverType = FEM_FULL; 
          if (glfemGetKey('B'))   newSolverType = FEM_BAND; 
          if (glfemGetKey('I'))   newSolverType = FEM_ITER; 
          if (glfemGetKey('X'))   newRenumType  = FEM_XNUM; 
          if (glfemGetKey('Y'))   newRenumType  = FEM_YNUM; 
          if (glfemGetKey('N'))   newRenumType  = FEM_NO; 

            
      
    } while(!glfemWindowShouldClose());
        
    glfemWindowFree();
    femDiffusionFree(theProblem);
    exit(EXIT_SUCCESS);
    return 0;
    
    

}

