#include "fem.h"


# ifndef NOEXPAND

void edgesExpand(femEdges *theEdges)
{
    // parcourt mesh pour obtenir les segments
    for(int i=0; i<theEdges->mesh->nElem; i++){
        for(int j=0; j<3; j++){
            // conditions pour pouvoir aller rechercher le noeud de départ
            if(j<2){
                (theEdges->edges[i*3+j]).node[0] = (theEdges->mesh)->elem[i*3+j];
                (theEdges->edges[i*3+j]).node[1] = (theEdges->mesh)->elem[i*3+j+1];
            }
            else{
                (theEdges->edges[i*3+j]).node[0] = (theEdges->mesh)->elem[i*3+j];
                (theEdges->edges[i*3+j]).node[1] = (theEdges->mesh)->elem[i*3];
            }
            // élément auquel appartient
            (theEdges->edges[i*3+j]).elem[0] = i;
            // valeur par défaut
            (theEdges->edges[i*3+j]).elem[1] = -1;
        }  
    }
}

# endif
# ifndef NOSORT

void edgesSort(femEdges *theEdges)
{
    qsort(theEdges->edges, theEdges->nEdge, sizeof(femEdge), edgesCompare);
}

# endif
# ifndef NOCOMPARE

# define MAX(a,b) a>b?a:b
# define MIN(a,b) a<b?a:b

int edgesCompare(const void* e0, const void *e1)
{   
    // getting min of e0
    int min1 = MIN(((femEdge*) e0)->node[0],((femEdge*) e0)->node[1]);

    // getting min of e1
    int min2 = MIN(((femEdge*) e1)->node[0],((femEdge*) e1)->node[1]);
    
    int diff = min1-min2;
    if(diff !=0) return diff;
    else{
        // getting max of e0
        int max1 = MAX(((femEdge*) e0)->node[0],((femEdge*) e0)->node[1]);

        // getting max of e1
        int max2 = MAX(((femEdge*) e1)->node[0],((femEdge*) e1)->node[1]);
        return max1-max2;
    }
}

# endif
# ifndef NOSHRINK

void edgesShrink(femEdges *theEdges)
{
    femEdge* edges_runner = theEdges->edges;
    int n = theEdges->nEdge;
    int nBoundary = theEdges->nBoundary;
    int nRunner = n;
    int i = 1;
    while(i<nRunner){
        edges_runner[i-nRunner+n] = edges_runner[i];
        int diag1 = edges_runner[i-1-nRunner+n].node[0]==edges_runner[i-nRunner+n].node[0];
        int diag2 = edges_runner[i-1-nRunner+n].node[1]==edges_runner[i-nRunner+n].node[1];
        int diag3 = edges_runner[i-1-nRunner+n].node[0]==edges_runner[i-nRunner+n].node[1];
        int diag4 = edges_runner[i-1-nRunner+n].node[1]==edges_runner[i-nRunner+n].node[0];
        if((diag1 & diag2) || (diag3 & diag4)){
            edges_runner[i-1-nRunner+n].elem[1] = edges_runner[i-nRunner+n].elem[0];
            n--;
            nBoundary-=2; 
        }
        i++;
    }
    
    // Reallocation du tableau des edges
    
    theEdges->edges = realloc(theEdges->edges, n * sizeof(femEdge));
    theEdges->nEdge = n;
    theEdges->nBoundary = nBoundary;
}

# endif
# ifndef NOBOUNDARYLENGTH

double edgesBoundaryLength(femEdges *theEdges)
{   
    double L = 0.0;
    for(int i=0; i<theEdges->nEdge; i++){
        if((theEdges->edges[i]).elem[1]==-1){
            // coordonnée du premier noeud
            double x1 = (theEdges->mesh)->X[theEdges->edges[i].node[0]];
            double y1 = (theEdges->mesh)->Y[theEdges->edges[i].node[0]];

            // coordonnée du deuxième noeud
            double x2 = (theEdges->mesh)->X[theEdges->edges[i].node[1]];
            double y2 = (theEdges->mesh)->Y[theEdges->edges[i].node[1]];

            // merci pythagore
            L += sqrt((x1-x2)*(x1-x2) +(y1-y2)*(y1-y2));
        }
    }
    return L;
}

# endif