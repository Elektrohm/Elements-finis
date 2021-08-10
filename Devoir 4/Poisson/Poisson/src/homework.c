#include"fem.h"
// Théo Denis 27411800


# ifndef NOPOISSONCREATE

femPoissonProblem *femPoissonCreate(const char *filename)
{
    femPoissonProblem *theProblem = malloc(sizeof(femPoissonProblem));
    theProblem->mesh  = femMeshRead(filename);           
    theProblem->edges = femEdgesCreate(theProblem->mesh);  
    if (theProblem->mesh->nLocalNode == 4) {
        theProblem->space = femDiscreteCreate(4,FEM_QUAD);
        theProblem->rule = femIntegrationCreate(4,FEM_QUAD); }
    else if (theProblem->mesh->nLocalNode == 3) {
        theProblem->space = femDiscreteCreate(3,FEM_TRIANGLE);
        theProblem->rule = femIntegrationCreate(3,FEM_TRIANGLE); }
    theProblem->system = femFullSystemCreate(theProblem->mesh->nNode);
    return theProblem;
}

# endif
# ifndef NOPOISSONFREE

void femPoissonFree(femPoissonProblem *theProblem)
{
    femFullSystemFree(theProblem->system);
    femIntegrationFree(theProblem->rule);
    femDiscreteFree(theProblem->space);
    femEdgesFree(theProblem->edges);
    femMeshFree(theProblem->mesh);
    free(theProblem);
}
    

# endif
# ifndef NOMESHLOCAL


void femMeshLocal(const femMesh *theMesh, const int i, int *map, double *x, double *y)
{
   for(int j=0; j<theMesh->nLocalNode; j++){
       map[j] = theMesh->elem[i*theMesh->nLocalNode+j];
       x[j] = theMesh->X[map[j]];
       y[j] = theMesh->Y[map[j]];
   }
}

# endif
# ifndef NOPOISSONSOLVE


void femPoissonSolve(femPoissonProblem *theProblem)
{
    int n = theProblem->space->n;

    // allocation des vecteurs utiles 
    double x[n], y[n], phi[n], dphidx[n], dphidy[n], dphidxsi[n], dphideta[n];
    int map[n];

    // boucle sur les i éléments du maillages
    for(int i=0; i<theProblem->mesh->nElem; i++){
        // récupère valeur élément local
        femMeshLocal(theProblem->mesh, i, map, x,y);

        // changement de variable
        for(int j=0; j<theProblem->rule->n; j++){
            // récupère paramètre intégration élément parents
            double xsi = theProblem->rule->xsi[j];
            double eta = theProblem->rule->eta[j];
            double weight = theProblem->rule->weight[j];

            // récupère fonction de forme
            femDiscretePhi2(theProblem->space, xsi, eta, phi);
            femDiscreteDphi2(theProblem->space, xsi, eta, dphidxsi, dphideta);


            // aide pour le changement de variable
            double dxdxsi = 0;
            double dxdeta = 0;
            double dydxsi = 0;
            double dydeta = 0;

            // variable changement de var 
            double jac = 0;

            for(int k=0; k<theProblem->space->n; k++){
                dxdxsi += x[k]*dphidxsi[k];
                dxdeta += x[k]*dphideta[k];

                dydxsi += y[k]*dphidxsi[k];
                dydeta += y[k]*dphideta[k];
            }

            jac = fabs(dxdxsi*dydeta-dydxsi*dxdeta);

            // applique le changement de variable
            for(int l=0; l<theProblem->space->n; l++){
                dphidx[l] = (1/jac)*(dydeta*dphidxsi[l]- dydxsi*dphideta[l]);
                dphidy[l] = (1/jac)*(dxdxsi*dphideta[l]- dxdeta*dphidxsi[l]);

            }
            // calcul élément [v,w] de A
            for(int v=0; v<theProblem->space->n; v++){
                for(int w=0; w<theProblem->rule->n; w++){
                    theProblem->system->A[map[v]][map[w]] += jac*weight*(dphidx[v]*dphidx[w] + dphidy[v]*dphidy[w]);
                }
            }
            // calcul élément [ib] de B
            for(int ib=0; ib<theProblem->space->n; ib++){
                theProblem->system->B[map[ib]] += phi[ib]*jac*weight;
            }
        }

    }

    for(int ic=0; ic<theProblem->edges->nEdge; ic++){
            if((theProblem->edges->edges[ic]).elem[1]==-1){
                femFullSystemConstrain(theProblem->system, (theProblem->edges->edges[ic]).node[0],0);
                femFullSystemConstrain(theProblem->system, (theProblem->edges->edges[ic]).node[1],0);
            }
        }
    // femFullSystemPrint(theProblem->system);
    theProblem->system->B = femFullSystemEliminate(theProblem->system);
}

# endif