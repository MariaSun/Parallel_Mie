//compile as : mpicc $(gsl-config --cflags) photon.c vector.c mesh.c mie_scat.c nrutil.c complex.c $(gsl-config --libs) -o photon

#include <mpi.h>
#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <time.h>

#include "mie_scat.h"
#include "complex.h"
#include "vector.h"
#include "nrutil.h"
#include "mesh.h"
#include </Users/Maria/gsl/include/gsl/gsl_sf_bessel.h>
#include </Users/Maria/gsl/include/gsl/gsl_complex.h>
#include </Users/Maria/gsl/include/gsl/gsl_complex_math.h>
#include </Users/Maria/gsl/include/gsl/gsl_sf_coupling.h>

#define	pi 3.1415926
#define a 0.01160239535139917 //mum  - scattering center radius



//////////////////////////////////////////////////////
////initialization of photon state////////////////////
//////////////////////////////////////////////////////

double theta_k=0.1;
double x=0.1;
int m_gamma=2, lambda=1;
double ux_p, uy_p, s_p, ell=0.;
double albedo=.9;
unsigned long nang=1; //remnance of Mie-scattering routine
double cxref_r=1.33, cxref_i=0.1;   //medium characteristic which never updates
double k_perp, k_z;

//////////////////////////////////////////////////////
////function declaration//////////////////////////////
//////////////////////////////////////////////////////

void distribution(double Theta, fcomplex *S1,  fcomplex *S2, double *S1_abs, double *S2_abs, double *p);
void find_max_p(double *p_max);
double photon(double p_max, double ux, double uy, double ballistic_length,double *X, double *Y, double *Int, double *weight);
void dist_mesh(double x[], double y[], double param[], int a_size, int mesh_size);

void field_gen(double ux, double uy, fcomplex *Ex, fcomplex *Ey);

int main(int argc, char **argv)
{
    int size, rank;
    MPI_Init(&argc, &argv);
    MPI_Status status;
    MPI_Comm_size(MPI_COMM_WORLD, &size);
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    
    int Nphotons, jj=0;
    double weight;
    
    //////////////////////////////////////////////////////
    ////read-off the beam flux////////////////////////////
    //////////////////////////////////////////////////////
    
    if( argc == 2 ) {
        Nphotons=atoi(argv[1]);
    }
    else if( argc > 2 ) {
        printf("Too many arguments supplied.\n");
        exit(EXIT_SUCCESS);
    }
    else if( argc < 2 ) {
        printf("Too fiew arguments supplied; 4 expected\n");
        exit(EXIT_SUCCESS);
    }
    else {
        printf("One argument expected.\n");
        exit(EXIT_SUCCESS);
    }
    
    if(Nphotons % (size-1) != 0){printf("number of photons %d should be divisible by number of threads used-1 %d\n", Nphotons, size-1); exit(1);}
    
    int step=Nphotons/(size-1);
    
    double Int;
    double p_max=0.;
    double ux, uy;
    double ballistic_length=1.;
    double ball_step=1.;
    k_perp=(x/a)*sin(theta_k), k_z=(x/a)*cos(theta_k);
    
    int output_size=5;
    
    
    for(int ii=1; ii<100; ii++){
        ballistic_length=ball_step*ii;
        if(rank==0)printf("%f ",ballistic_length);
    
    double *bufs = NULL;
    double *buf = NULL;
    Vector pre_buf;
    vector_init(&pre_buf);
    int *size_survived=NULL, count=0, out_size=(output_size * Nphotons)+output_size,*displs = NULL,full_size=0;
    
    /////////////////////////////////////////
    ////generating angles////////////////////
    /////////////////////////////////////////
    
    if(rank==0)
    {
        count=output_size;
        bufs=malloc(sizeof(double) * out_size);
        
        buf=malloc(sizeof(double) * output_size);
        for(int ii=0;ii<output_size; ii++) buf[ii]=0.;
        
        find_max_p(&p_max);
        double count=1;
        do{
            MPI_Send(&p_max, 1, MPI_DOUBLE, size-count, 1, MPI_COMM_WORLD);
            count++;}while(count<size);
    }
    else{MPI_Recv(&p_max, 1, MPI_DOUBLE, 0, 1, MPI_COMM_WORLD, &status);}//communicate max of scattering distribution for the current parameters of the scattering centers
    MPI_Barrier(MPI_COMM_WORLD);
    
    
    //Form the vector buf[..] that contains all the information about the photons that are coming though
    if(rank!=0){
    count=0;
    srand48((unsigned) time(NULL)*rank);//rank
    for(jj=0; jj<step; jj++)
        {
            //generating a circular incoming beam
            double phi=2*pi*drand48(), r;
            double u=drand48()+drand48();
            if(u>1)r=2-u; else r=u;
            double ux_init=r*cos(phi), uy_init=r*sin(phi);
            
            photon(p_max, ux_init, uy_init, ballistic_length,&ux, &uy, &Int, &weight);
            if(weight>1.E-8){
                vector_set(&pre_buf, output_size*count, ux);
                vector_set(&pre_buf, (output_size*count)+1, uy);
                vector_set(&pre_buf, (output_size*count)+2, Int);
                vector_set(&pre_buf, (output_size*count)+3, weight);
                
                double dist=pow(ux*ux+uy*uy,.5);//distance from the center of the beam
                vector_set(&pre_buf, (output_size*count)+4, dist);
                count++;
            }
            
        }
        
        if(count>0){

            buf=malloc(sizeof(double) * (count*output_size));

            for(int ii=0;ii<count;ii++){
                for(int jj=0;jj<output_size;jj++){
                    int cc=(output_size*ii)+jj;
                    buf[cc]=vector_get(&pre_buf, cc);
                }
            }
        }else buf=malloc(sizeof(double));

        count=count*output_size;
    }
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //prepare the displs[..] - displacement vector for MPI_Gatherv
    if (rank==0){
        displs = malloc( (size+1) * sizeof(int) );
        size_survived = malloc(size * sizeof(int));
        for(int ii=0; ii<size; ii++){
            displs[ii]=0;
        }
    }
    
    //collects the number of survived photons times the number of parameters per photon from all threads
    //assigns the displacement vector
    MPI_Gather(&count, 1, MPI_INT, size_survived, 1, MPI_INT, 0, MPI_COMM_WORLD);
    if(rank==0){
        full_size=0;
        displs[0]=0;
        int tmp[size];
        for(int jj=0; jj<size; jj++)tmp[jj]=size_survived[jj];
        for(int ii=1;ii<=size; ii++){
            for(int jj=0; jj<ii; jj++)displs[ii]+=tmp[jj];
        }
        
        for(int ii=0; ii<size; ii++){
            full_size+=tmp[ii];
        }
    }

    MPI_Barrier(MPI_COMM_WORLD);
    
    //collect all the info about the survived photon onto master thread
    MPI_Gatherv(buf, count, MPI_DOUBLE, bufs, size_survived, displs, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    
    MPI_Barrier(MPI_COMM_WORLD);
    
    //the information about the photons is on the master thread now
    if(rank==0)
    {
        double *photon_info[6], *final_set[6];
        int photon_number=full_size/output_size-1;
    
        for(int ii=0; ii<6; ii++)photon_info[ii]=malloc(sizeof(double)*photon_number);
        for(int ii=0; ii<6; ii++)final_set[ii]=malloc(sizeof(double)*photon_number);
        
        for(int ii=0; ii<photon_number; ii++)photon_info[0][ii]=ii;//index
        for(int ii=output_size; ii<full_size; ii++){
            int index=(int)(ceil(ii/output_size)-1);
            if(ii%output_size==0) {photon_info[2][index]=bufs[ii];}
            if(ii%output_size==1) {photon_info[3][index]=bufs[ii];}
            if(ii%output_size==2) photon_info[4][index]=bufs[ii];
            if(ii%output_size==3) photon_info[5][index]=bufs[ii];
            if(ii%output_size==4) photon_info[1][index]=bufs[ii];
        }
        
        //sorts out top 10 photons by distance and throughs those away
        int top_percent=(full_size/100)*10;
        for(int ii=0; ii<photon_number; ii++){
            
            for(int jj=0; jj<photon_number; jj++){
                if(photon_info[1][jj]>photon_info[1][ii]){
                    double tmp[5];
                    for(int kk=0; kk<6; kk++){
                        tmp[kk] = photon_info[kk][ii];
                        photon_info[kk][ii] = photon_info[kk][jj];
                        photon_info[kk][jj]=tmp[kk];
                    }
                }
            }
        }
        
        
        //forms final set of photons, which are going to be meshed and analysed
        int newphoton_number=photon_number-top_percent;
        
        for(int ii=0; ii<newphoton_number; ii++){
            for(int jj=0; jj<6; jj++){
                final_set[jj][ii]=photon_info[jj][ii];
            }
        }
        
        dist_mesh(photon_info[2], photon_info[3], photon_info[5], newphoton_number, 100);
        
        
//        FILE *row_data;
//        row_data=fopen("test.txt", "w");
//        for(int ii=output_size; ii<full_size; ii++)
//        {
//            double current_value=bufs[ii];
//                fprintf(row_data,"%.10f ",bufs[ii]);
//        };
//        fclose(row_data);
        
        for(int ii=0; ii<6; ii++)free(photon_info[ii]);
        for(int ii=0; ii<6; ii++)free(final_set[ii]);
        
    }
    
    
    free(buf);
    
    if(rank!=0)vector_free(&pre_buf);
    if(rank==0){
        free(bufs);
        free(size_survived);
        free(displs);
    }
    }
    MPI_Finalize();
    return 0;

}

double photon(double p_max, double ux, double uy, double ballistic_length,double *X, double *Y, double *Int, double *weight)
/////////////////////////////////////////////////
////photon MC simulation/////////////////////////
////suppy: max of scattering distribution////////
/////////////////////////////////////////////////
{
    double Phi, Theta, S1_abs, S2_abs, p=0.,f, w, F, s=0.;
    fcomplex S1=Complex(0.,0.), S2=Complex(0.,0.);
    fcomplex Ex_p, Ey_p;
    
        w=1.;
    
        
        fcomplex Ex=Complex(1,0); //(Ex==E_parallel and Ey==E_perp)
        
        fcomplex Ey=Complex(0,0);

    
        do
        {
        negativez:
            
            do
            {
                Theta = pi*drand48()*2;
                
                f=drand48()*p_max;
                distribution(Theta, &S1, &S2, &S1_abs, &S2_abs, &p);
            }
            while(f>p);
            
            Phi = 2*drand48()*pi;
            
            /////////////////////////////////////////
            ////update photon weight/////////////////
            /////////////////////////////////////////
            
            w=w*albedo;//I think it is "-" instead
            
            if (w<=0.00000000001) {
                *weight=w;
                
                *X=0.;
                *Y=0.;
                
                *Int=0.;
                goto exit;}//photon is dead
            
            /////////////////////////////////////////
            ////update position//////////////////////
            /////////////////////////////////////////
            
            
            ux_p = cos(Theta)*cos(Phi)*ux + cos(Theta)*sin(Phi)*uy - sin (Theta)*s;
            uy_p = -sin(Phi)*ux + cos(Phi)*uy;
            s_p = sin(Theta)*cos(Phi)*ux + sin(Theta)*sin(Phi)*uy + cos (Theta)*s;
            
            
            ell+=pow((ux_p-ux)*(ux_p-ux)+(uy_p-uy)*(uy_p-uy)+(s_p-s)*(s_p-s),.5);//optical path
            
            ux+=ux_p;
            uy+=uy_p;
            s+=s_p;
            
            if(s<=0.) goto negativez;//unphysical step
            
            /////////////////////////////////////////
            ////update fields////////////////////////
            /////////////////////////////////////////
            
            distribution(Theta, &S1, &S2, &S1_abs, &S2_abs, &p);
            
            fcomplex tt=Cmul(Ex,Conjg(Ey));
            
            F=pow(((S2_abs*S2_abs*cos(Phi)*cos(Phi)+S1_abs*S1_abs*sin(Phi)*sin(Phi))*Cabs(Ex)*Cabs(Ex) + (S2_abs*S2_abs*sin(Phi)*sin(Phi)+S1_abs*S1_abs*cos(Phi)*cos(Phi))*Cabs(Ey)*Cabs(Ey) + 2*(S2_abs*S2_abs - S1_abs*S1_abs)*cos(Phi)*sin (Phi)*tt.r),-.5);
            
            fcomplex B11, B12, B21, B22;
            
            B11=Cmul(S2,Complex(cos(Phi),0.));
            B12=Cmul(S2,Complex(sin(Phi),0.));
            B21=Cmul(S1,Complex(sin(Phi),0.));
            B22=Cmul(S1,Complex(cos(Phi),0.));
            
            Ex_p = RCmul(F,Cadd(Cmul(B11,Ex), Cmul(B12,Ey)));
            Ey_p = RCmul(F,Csub(Cmul(B22,Ey),Cmul(B21,Ex)));
            
            if(Ex.r==0&&Ex.i==0&&Ey.r==0&&Ey.i==0) goto exit;//photon is dead
            
            Ex=Cmul(Complex(Ex_p.r,Ex_p.i),Complex(ux,0.));
            Ey=Cmul(Complex(Ey_p.r,Ey_p.i),Complex(uy,0.));
            
        }while(s<=ballistic_length);//thickness
        
        Ex=Cmul(Ex,Complex(cos((x/a)*ell),sin((x/a)*ell)));
        Ey=Cmul(Ey,Complex(cos((x/a)*ell),sin((x/a)*ell)));
    
    *X=ux;
    *Y=uy;

    *Int=pow(w,-.5)*(double)(Cabs(Ex)*Cabs(Ex)+Cabs(Ey)*Cabs(Ey));
    *weight=w;
    
    exit:
    
    
    return 0;
    
}

void find_max_p(double *p_max)
/////////////////////////////////////////////////
////generates max of p(theta)////////////////////
////suppy: scatterers parameters/////////////////
/////////////////////////////////////////////////
{
    
    //this function needs to be tested, the search algorythm is very primitive
    
    
    int i,j;
    double p_temp1=0.,p_temp2=0., theta_min=0., theta_max=pi, Theta=0., p=0., delta=0., S1_abs, S2_abs;
    fcomplex S1, S2;
    
    do
    {
        Theta=(theta_max+theta_min)/2;
        
        for(i=1;i<=100;i++){
            
            delta=(Theta-theta_min)/100;
            distribution(theta_min+(i*delta), &S1, &S2, &S1_abs, &S2_abs, &p);
            p_temp1+=p;
        }
        p_temp1=p_temp1/100;
        
        for(j=1;j<=100;j++){
            
            delta=(theta_max-Theta)/100;
            distribution(Theta+(i*delta), &S1, &S2, &S1_abs, &S2_abs, &p);
            p_temp2+=p;
        }
        p_temp2=p_temp2/100;
        
        if(p_temp1>p_temp2) theta_max=Theta;
        else theta_min=Theta;
       
    }while (fabs(p_temp1-p_temp2)>0.00001);
    
    if(p_temp1>p_temp2) *p_max=p_temp1;
    else *p_max=p_temp2;
    
    
}


void distribution(double Theta, fcomplex *S1,  fcomplex *S2, double *S1_abs, double *S2_abs, double *p)
    ///////////////////////////////////////////////////////
    ////generates distribution p(theta)////////////////////
    ////suppy: scatterers parameters, theta////////////////
    ///////////////////////////////////////////////////////
{
    double cxs1_r, cxs1_i, cxs2_r, cxs2_i;
    float qext, qsca, qback, gsca;
    
    void bhmie(float x, double cxref_r, double cxref_i, unsigned long nang, double Theta, double *cxs1_r, double *cxs1_i, double *cxs2_r, double *cxs2_i, float  *qext, float *qsca, float *qback, float *gsca);
    
    bhmie(x, cxref_r, cxref_i, nang, Theta, &cxs1_r, &cxs1_i, &cxs2_r, &cxs2_i, &qext, &qsca, &qback, &gsca);
    
    
    *S1_abs = cxs1_r*cxs1_r + cxs1_i*cxs1_i;
    *S2_abs = cxs2_r*cxs2_r + cxs2_i*cxs2_i;
    
    *S1 = Complex(cxs1_r, cxs1_i);
    *S2 = Complex(cxs2_r, cxs2_i);
    
    *p=(cxs1_r*cxs1_r + cxs1_i*cxs1_i + cxs2_r*cxs2_r + cxs2_i*cxs2_i)/(qsca*x*x);
}
