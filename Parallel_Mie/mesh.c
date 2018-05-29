//M.Solyanik 05/28/2018///////////////////////////////////////////////////////
//The function does meshing by mesh_size - steps of the supplied x[] and y[]
//by addition of the parameter param[]
//a_size - size of arrays x[] and y[]

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>

void dist_mesh(double x[], double y[], double param[], int a_size, int mesh_size){
    
    double max_x=0., min_x=0., max_y=0., min_y=0., dmesh_x, dmesh_y;
    int N_x, N_y;
    double pre_mesh[3][mesh_size][mesh_size];

    
    for(int i=0; i<mesh_size; i++){
        for(int j=0; j<mesh_size; j++){
            for(int k=0; k<3; k++)pre_mesh[k][i][j]=0.;
        }
    }
    
    for(int i=0; i<a_size; i++){
        if(x[i]>max_x) max_x=x[i];
        if(x[i]<min_x) min_x=x[i];
        
        if(y[i]>max_y) max_y=y[i];
        if(y[i]<min_y) min_y=y[i];
    }
    
    dmesh_x=fabs((max_x-min_x)/mesh_size);
    dmesh_y=fabs((max_y-min_y)/mesh_size);
    
    for(int i=0; i<a_size; i++){
        
        if(x[i]<0){
            N_x=(int)fabs(ceil((min_x-x[i])/dmesh_x));
        }else{
            N_x=(int)floor(x[i]/dmesh_x) +fabs(min_x/dmesh_x);
        }
        
        if(y[i]<0){
            N_y=(int)fabs(ceil((min_y-y[i])/dmesh_y));
        }else{
            N_y=(int)floor(y[i]/dmesh_y) +fabs(min_y/dmesh_y);
        }
        
        pre_mesh[2][N_x][N_y]+=param[i];
        pre_mesh[0][N_x][N_y]=min_x + (dmesh_x*(N_x+.5));
        pre_mesh[1][N_x][N_y]=min_y + (dmesh_y*(N_y+.5));
        
    }
    
    printf("%f ", (max_x+max_y)/2);
    
}

