#include "glfem.h"

int main(void)
{  
    femMesh  *theMesh  = femMeshRead("../data/gear60.txt");    
    femEdges *theEdges = femEdgesCreate(theMesh);    
    
    edgesExpand(theEdges);               //   femEdgesPrint(theEdges);
    edgesSort(theEdges);                 //   femEdgesPrint(theEdges);
    edgesShrink(theEdges);               //   femEdgesPrint(theEdges);
    printf("Boundary edges  : %i \n", theMesh->nNode);
    printf("Boundary length : %14.7e \n", edgesBoundaryLength(theEdges));

    char theMessage[256];
    sprintf(theMessage, "Boundary edges : %i", theEdges->nBoundary);
      
//
//  On superpose le maillage (en bleu), 
//  tous les segments frontieres (en noir),
//  et la frontiere (en rouge)
//
//  Au depart de votre travail, vous devriez obtenir un maillage bleu....
//  et a la fin de l'exercice un maillage noir avec bord rouge :-)
//

    glfemWindowCreate("EPL1110 : Edges",480,480,theMesh->nNode,theMesh->X,theMesh->Y);
    do
    {
        glfemWindowUpdate();    
        
        glfemSetLineWidth(0.0001); 
        glfemSetColor(GLFEM_BLUE); 	  glfemPlotMesh(theMesh);  
        glfemSetLineWidth(0.0005); 
        glfemSetColor(GLFEM_BLACK); 	glfemPlotEdges(theEdges);  
        glfemSetLineWidth(0.0015);
        glfemSetColor(GLFEM_RED); 	  glfemPlotBnd(theEdges);  
       
        glfemDrawMessage(theMessage,(double[2]){16.0, 30.0}); 
    } while(!glfemWindowShouldClose());

    glfemWindowFree();
    femEdgesFree(theEdges);
    femMeshFree(theMesh);

    exit(EXIT_SUCCESS);
    return 0;
    
 
}



