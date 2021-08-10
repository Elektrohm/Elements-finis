
#include "glfem.h"

double inertiaGearIntegrate(femMesh *theMesh, femIntegration *theRule);

int main(void)
{  
    femIntegration* theRule = femIntegrationCreate(3, FEM_TRIANGLE);
    femMesh *theMesh = inertiaGearMeshRead("../data/gear60.txt");   
           
    double rho = inertiaGearSteelRho();
    double   I = inertiaGearInertia(theMesh,theRule,rho); 

    double *theField = malloc(sizeof(double)*theMesh->nNode); int i;
    for (i=0; i < theMesh->nNode; ++i) {
        double xLoc = theMesh->X[i];
        double yLoc = theMesh->Y[i];
        theField[i]  = xLoc*xLoc + yLoc*yLoc; }
           
    char theMessage[256];       
    sprintf(theMessage,"Inertia = %14.7e kg m^2",I);
    printf("%s\n",theMessage);

    glfemWindowCreate("EPL1110 : Inertia of Gear",480,480,theMesh->nNode,theMesh->X,theMesh->Y);
    do
    {
        glfemWindowUpdate(); 
        glfemSetColor(GLFEM_BLACK);  
        glfemPlotSolution(theMesh,theField);    
        
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0}); 
    } while(!glfemWindowShouldClose());

    free(theField);
    femIntegrationFree(theRule);
    inertiaGearMeshFree(theMesh);
    exit(EXIT_SUCCESS);    
    return 0; 
 
}





 
