#include <stdio.h>
#include <math.h>

#include "fem.h"
#include "glfem.h"



int main(void)
{   
 
    femPoissonProblem* theProblem = femPoissonCreate("../data/triangles22.txt");
    
    // Pour Windows, remplacer l'argument :
    // ("../data/triangles_166.txt") 
    // par :
    // ("..\\data\\triangles_166.txt") 
    //
    // Sorry for the inconvenience :-)
    // On réfléchit pour rendre cela plus transparent dans les homeworks suivants :-)
    // Be patient !
     
    printf("Number of elements    : %4d\n", theProblem->mesh->nElem);
    printf("Number of local nodes : %4d\n", theProblem->mesh->nLocalNode);
    printf("Number of segments    : %4d\n", theProblem->edges->nBoundary);
    printf("Number of unknowns    : %4d\n", theProblem->system->size);

    femPoissonSolve(theProblem);  
  
 
    printf("Maximum value : %.4f\n", femMax(theProblem->system->B,theProblem->system->size));
    fflush(stdout);
    
    char theMessage[256];
    sprintf(theMessage, "Max : %.4f", femMax(theProblem->system->B,theProblem->system->size));
    
    glfemWindowCreate("EPL1110 : Poisson",480,480,theProblem->mesh->nNode,theProblem->mesh->X,theProblem->mesh->Y);
    do {
      glfemWindowUpdate();    
        
      glfemPlotSolution(theProblem->mesh,theProblem->system->B);
      glfemPlotMesh(theProblem->mesh);  

    
      glfemDrawMessage(theMessage,(double[2]){16.0, 30.0});    
    } while(!glfemWindowShouldClose());
        
    femPoissonFree(theProblem);
    glfemWindowFree();
    exit(EXIT_SUCCESS);
    return 0;
    
  

}

