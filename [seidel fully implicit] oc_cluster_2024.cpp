#include <cstring>
#include <stdio.h>
#include <iostream>
#include <string.h>
#include <cstring>
#include <math.h>
#include <fstream>
#include <time.h>
#include <algorithm>

/*
OVARIAN CANCER CELL DYNAMICS:
[!!!] RUNNING THE CODE NEEDS SUBFOLDER "numerical_results"
[!!!] THIS CODE BELOW RUNS 20 REALIZATIONS
*/


const int Number_of_Realizations = 20;
const double L=20.0;// linear size of the substratea
const int M=101; // grid size is M x M
const int Nc=32; // number of cells
const double dt=0.005;//0.0001; // time step
const int Nt=10000;//10000;//100000; //70000; // number of time steps
const double E0=0.5; // elastic coefficient for rho=0
const double E1=2.5; // elastic coefficient at rho
const double gamma_1=0.01;//0.1;//<-This one is correct//0.01;//0.025; // rate of degradation due to cells
const double gamma_2=0.002;//0.02//<-This one is correct;
const double beta_0=0.4;//0.0125; // rate of self-degradation of chemical
const double v_prop=1.0;//0.1;//0.2;//1.0;//0.03;//0.05; //0.1// crawling velocity
const double k_rho=0.0; // proportionality coefficient between velocity and gradient of surface
const double C_0=1.0;//0.00125; //rate of secretion
const double d_cell=1.0;// diameter of cell
const double k_int=0.5; // <----------------------------------------------------ZEROED  interaction coefficient
const double k_sigma=-1.0; //20.0 reorientation towards principal axis of stress
const double D_C=0.1;//1.0;//0.25; // spatial diffusion
const double D_phi=0.5;//0.15; // rotational diffusion
const double eta=0.125;//3.0; //viscosity
const double ell_shoulder=L/(M*1.0)*sqrt(5.0); //cell shoulder length (sqrt(5))
const double alpha_diss=4.0;
const double weak_adhesion=0; // 0.05// adhesion between cells of type 2
const double normal_adhesion=0;// 100.0 adhesion between cells of type 1
const double intermediate_adhesion=0.0;// 1.0 adhesion between cells of different types
const double distance_adhesion=2.25;
//const double k_u=3.0; //how substrate displacement affect cells
const double zeta=0.5; //0.25 traction coefficient
const int record_frequency=100;//20;//250;


double u1[M][M],u1_old[M][M],u2[M][M],u2_old[M][M],u1_iter_old[M][M],u2_iter_old[M][M];
double b1_seidel[M*M],b2_seidel[M*M];//substrate deformations
double load[M][M],load_old[M][M];// load on substrate
double stress_xx[M][M],stress_xy[M][M],stress_yy[M][M],magnitude_of_displacement[M][M];// substrate stress
double dlambda[M][M];
double phi_stress[M][M]; //stress principal orientation angle
double traction_x[M][M],traction_y[M][M];// traction field
double rho[M][M],rho_old[M][M];//density of cross-links
double chemical[M][M],chemical_old[M][M],chemical_iter_old[M][M];//cell chemical distribution
double b_chemical[M*M];//right-hand side for seidel chemical
double x[Nc],x_old[Nc];// x-coordinates of cells
double y[Nc],y_old[Nc];// y-coordinates of cells
double v1[Nc],v2[Nc];
double phi[Nc],phi_old[Nc];// phi-coordinate of cells
double u_right[Nc],u_left[Nc];// side propulsion to compute torque


//continuous periodicity
double per(double x)
{
    double helper=x;
    if (x>L) helper=helper-L;
    if (x<0) helper=helper+L;
    return helper;
}

//continuous periodicity - recursivity idea
double per2(double x)
{
   if ((-L/2<x)&&(x<=L/2)) return x;
   if (x>L/2) return per2(x-L);
   if (x<=-L/2) return per2(x+L);
}

//continuous periodicity in [0;L] - no recursive
double per3(double x)
{
    int helper_int;
    if (x==0) return 0.0;
    if (x>0) helper_int=(int)(x/L);
    if (x<0) helper_int=(int)(x/L)-1;
    return x-(helper_int)*L;
}

//continuous periodicity in [-L/2;L/2] - no recursive
double per4(double x)
{
    int int_helper=(x>=0)-(x<0);
   return int_helper*(per3(fabs(x)+L/2.0)-L/2.0);
}


//index periodicity
int per_index(int i)
{
    int helper=i;
    if (i>M-1) helper=helper-M;
    if (i<0) helper=helper+M;
    return helper;
}


//bilinear interpolation
double interpol(double x,double y,double x_1,double y_1,double f_11,double f_12,double f_21,double f_22)
{

if (x<x_1) printf("Something is wrong %f < %f\n",x,x_1);

double dx=L/(1.0*M);
double X1=(x_1+dx-x)/dx;
double X2=(x-x_1)/dx;
double Y1=(y_1+dx-y)/dx;
double Y2=(y-y_1)/dx;

return f_11*X1*Y1+f_12*X2*Y1+f_21*X1*Y2+f_22*X2*Y2;
}


//young modulus
double young_modulus(double rho_1)//"1": one
{
    return E0+E1*rho_1;
}





double mass_of_cell(double x,double y,double x_cell,double y_cell)
{
    double r;// distance between cell and the point
    double l_0=1.0;// decay rate
    double F0=5.0;// amplitude
    r=sqrt(pow(per4(x_cell-x),2)+pow(per4(y_cell-y),2));

    return (r<2.0)*F0*exp(-pow(r/l_0,2));
}

double mass_of_cell_traction(double x,double y,double x_cell,double y_cell)
{
    double r;// distance between cell and the point
    double l_0=2.5;// decay rate
    double F0=2.5;//F0=5.0;// amplitude
    r=sqrt(pow(per4(x_cell-x),2)+pow(per4(y_cell-y),2));

    return (r<5.0)*F0*exp(-pow(r/l_0,2));
}

double cell_interaction(double x1,double y1,double x2,double y2,int k)//k=1,2: coordinate
{
    double r;// distance between two cells
    double helper;
    double x_per,y_per;
    x_per=per4(x1-x2);
    y_per=per4(y1-y2);
    r = sqrt(pow(x_per,2)+pow(y_per,2)); if (r==0) return 0.0;
    double F1=5.0;// 5.0//interaction strength
    double l_0=1.0;// 1.0
   // double l1=0.2;// radius of attraction
    //double l2=0.1;// radius of steric repulsion
    helper=(r<1.0)*F1*exp(-pow(r/l_0,2))/r;//r;
    if (k==1) return helper*x_per;
    if (k==2) return helper*y_per;
}

//standard max and min functions for integer variables
int min_int(int k_1, int k_2)
{
    if (k_1<k_2) return k_1; else return k_2;
}
int max_int(int k_1, int k_2)
{
    if (k_1<k_2) return k_2; else return k_1;
}

double min_double(double x_1,double x_2)
{
    if (x_1<x_2) return x_1; else return x_2;
}

int main()
{
//To work with files
char ch[80];
FILE *substrate;
FILE *cells;
FILE *loads;
FILE *chemical_file;
FILE *displacement_magnitude;


//initialization
int i,j,it,n,k; // indexes for cycles
double ran; // variable to be used for randomizations
double xi,yj; //for site's coordinates
double dx=L/(1.0*M); //spatial step - resolution of grid
double dload;
double x_interaction_nk,y_interaction_nk,scalar_interaction;
//double stress_xx,stress_yy,stress_xy; //rhs equation for u
double w_average;//average w (mass of cell) for computation of traction force
double E_local; // local young modulus
double discriminant; // parameters needed for computing for principal direction
double lambda_1,lambda_2; // see above
double sigma_1,sigma_2; // parameters to move from low stress to large
double gradient_abs_sigma_x,gradient_abs_sigma_y; // gradient of absolute value of sigma
double x_shoulder,y_shoulder;
double f_11,f_12,f_21,f_22;
int seidel_iteration,seidel_c_iteration;
double error_norm,error_c_norm;
double diag,u1ip1j,u1im1j,u1ijp1,u1ijm1,u1ip1jp1,u1im1jm1,u1ip1jm1,u1im1jp1,u2ip1j,u2im1j,u2ijp1,u2ijm1,u2ip1jp1,u2im1jm1,u2ip1jm1,u2im1jp1;
double diag_c,cip1j,cim1j,cijp1,cijm1,cip1jp1,cim1jm1,cip1jm1,cim1jp1;
double varkappa;
double max_u_kmkp1,max_u,helper_seidel,helper_seidel_2;
double max_c_kmkp1,max_c;
double force_adhesion_x,force_adhesion_y;
double helper_adhesion;
double r_distance;
int cell_label[Nc];
double r_x_per,r_y_per;
int R;
double ogamma_2;

for (R=1;R<Number_of_Realizations+1;R++)
{

ogamma_2 = gamma_2;

// open files
sprintf(ch,"numerical_results/substrate_density%d.txt",R);
substrate=fopen(ch,"w");
sprintf(ch,"numerical_results/cells%d.txt",R);
cells=fopen(ch,"w");
sprintf(ch,"numerical_results/loads%d.txt",R);
loads=fopen(ch,"w");
sprintf(ch,"numerical_results/chemical%d.txt",R);
chemical_file=fopen(ch,"w");
sprintf(ch,"numerical_results/displ%d.txt",R);
displacement_magnitude=fopen(ch,"w");

//computation of w_average
w_average=0.0;
for (i=0;i<M;i++)
{
    xi=(i+0.5)*dx;
    for (j=0;j<M;j++)
    {
        yj=(j+0.5)*dx;
        w_average=w_average+dx*dx*(mass_of_cell_traction(0.5*L,0.5*L,xi,yj));
    }
}

// CHECK IF THE CODE WORKS
printf("Am I still working?");
//labeling cells <--- Correct if you want two species or more than 100 cells!!!
for (i=0;i<Nc;i++)
{
    cell_label[i]=1;
}


printf("The Nc is %d",Nc);


//checking how per3 works
printf("per3(%f)=%f \n",-10.0,per3(-10.0));
printf("per3(%f)=%f  \n",-110.0,per3(-110.0));
printf("per3(%f)=%f \n",10.0,per3(10.0));
printf("per3(%f)=%f  \n",110.0,per3(110.0));
printf("per3(%f)=%f \n",310.0,per3(310.0));
printf("per3(%f)=%f  \n",-310.0,per3(-310.0));


//checking how per4 works
printf("per4(%f)=%f \n",-10.0,per4(-10.0));
printf("per4(%f)=%f  \n",-60.0,per4(-60.0));
printf("per4(%f)=%f \n",10.0,per4(10.0));
printf("per4(%f)=%f  \n",110.0,per4(110.0));
printf("per4(%f)=%f \n",310.0,per4(310.0));
printf("per4(%f)=%f  \n",-310.0,per4(-310.0));

srand(time(NULL));


// INITIALIZATION
for (n=0;n<Nc;n++)
{
    ran = ((double) rand() / (RAND_MAX));
    x[n]=0.5*L+L/2.0*(2.0*ran-1.0); //<-- TYPICAL
    //x[n]=0.5*L;
    /*if (n<3) x[n]=((n+0.5)/(3*1.0))*L;
    if ((n<6)&&(n>2)) x[n]=((n-3+0.5)/(3*1.0))*L;
    if (n>5) x[n]=((n-6+0.5)/(3*1.0))*L;*/

    ran = ((double) rand() / (RAND_MAX));

    /*if (n<3) y[n]=0.5*L;
    if ((n<6)&&(n>2)) y[n]=0.75*L;
    if (n>5) y[n]=0.25*L;*/

    y[n]=0.5*L+L/2.0*(2.0*ran-1.0); //<-- TYPICAL
    //y[n]=0.5*L-2.0;
    //0.5*L+L/2.0*(2.0*ran-1.0);
    //y[n]=L/2.0;
    ran = ((double) rand() / (RAND_MAX));
    phi[n]=2.0*M_PI*ran; //<-- TYPICAL
    //phi[n]=M_PI/2.0;
}

for (i=0;i<M;i++)
{
    for (j=0;j<M;j++)
    {
        load[i][j]=0.0;
        xi=(i)*dx;
        yj=(j)*dx;
        for (n=0;n<Nc;n++)
            {
                load[i][j]=load[i][j]+mass_of_cell(xi,yj,x[n],y[n]);
            }
        u1[i][j]=0.0;
        u2[i][j]=0.0;
        //rho[i][j]=3.0*sin(M_PI*yj/L)+0.25;//instead of 3.0 I had 0.75
        //rho[i][j]=-3.0*sin(M_PI*yj/L)+3.25;
        rho[i][j]=1.0;
        //rho[i][j] = (1 + sin(2*M_PI*xi/L))/2.0;
        chemical[i][j]=0.0;
    }
}

for (n=0;n<Nc;n++)
{
    printf("%f %f %f\n",mass_of_cell(dx/2,dx/2,x[n],y[n]),x[n],y[n]);
}
printf("dx=%f\n",dx);

//dynamics
for (it=0;it<Nt;it++)
{
    //new -> old
    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            u1_old[i][j]=u1[i][j];
            u2_old[i][j]=u2[i][j];
            load_old[i][j]=load[i][j];
            rho_old[i][j]=rho[i][j];
            chemical_old[i][j]=chemical[i][j];
        }
    }

    for (n=0;n<Nc;n++)
    {
        x_old[n]=per3(x[n]);
        y_old[n]=per3(y[n]);
        phi_old[n]=phi[n];
    }

    //computation of traction
    for (i=0;i<M;i++)
    {
        //xi=(i+0.5)*dx;
        xi=i*dx;
        for (j=0;j<M;j++)
        {
            //traction of cell
            yj=j*dx;

            traction_x[i][j]=0.0;
            traction_y[i][j]=0.0;
            for (n=0;n<Nc;n++)
            {
                traction_x[i][j]=traction_x[i][j]+zeta*(mass_of_cell_traction(xi,yj,x_old[n],y_old[n])-(mass_of_cell_traction(xi,yj,x_old[n],y_old[n])>0.001)*w_average)*cos(phi_old[n]);
                traction_y[i][j]=traction_y[i][j]+zeta*(mass_of_cell_traction(xi,yj,x_old[n],y_old[n])-(mass_of_cell_traction(xi,yj,x_old[n],y_old[n])>0.001)*w_average)*sin(phi_old[n]);
            }
        }
    }
    // Dynamics

    // Updating substrate properties: displacement,density,chemical;
    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            xi=(i+0.5)*dx;
            yj=(j+0.5)*dx;
            load[i][j]=0.0;

            for (n=0;n<Nc;n++)
            {
                // HERE DIFFERENCE IN SECRETION IS TAKEN INTO ACCOUNT
                load[i][j]=load[i][j]+(10.0*(cell_label[n]==1)+1.0*(cell_label[n]==2))*mass_of_cell(xi,yj,x_old[n],y_old[n]);
            }


        }
    }

    // Here we will update the displacement by using Gauss-Seidel algorithm
    // The system will be rewritten as Ax=b
    // Compute b
    varkappa=pow(dx,2)/dt;
    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            b1_seidel[j*M+i]=eta*u1_old[i][j]*varkappa+traction_x[i][j]*pow(dx,2);
            b2_seidel[j*M+i]=eta*u2_old[i][j]*varkappa+traction_y[i][j]*pow(dx,2);
        }
    }


    //iteration cycle
    seidel_iteration=0;
    error_norm=0.5;


    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            //u1[i][j]=0.0;
            //u2[i][j]=0.0;
            u1_iter_old[i][j]=1.0;
            u2_iter_old[i][j]=1.0;
        }
    }

    while ((seidel_iteration<1000)&&(error_norm>0.0001))
    {
        //count number of iterations
        seidel_iteration++;

        //count error_norm
        error_norm=0.0;

        max_u_kmkp1=0.0;
        max_u=0.0;
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                helper_seidel=sqrt(pow(u1[i][j]-u1_iter_old[i][j],2)+pow(u2[i][j]-u2_iter_old[i][j],2));
                helper_seidel_2=sqrt(pow(u1[i][j],2)+pow(u2[i][j],2));
                if (helper_seidel>max_u_kmkp1) max_u_kmkp1=helper_seidel;
                if (helper_seidel_2>max_u) max_u=helper_seidel_2;
            }
        }
        error_norm=helper_seidel/helper_seidel_2;
       // error_norm=1.0*(seidel_iteration<20);
        //old -> to -> new
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                u1_iter_old[i][j]=u1[i][j];
                u2_iter_old[i][j]=u2[i][j];
            }
        }

        //Gauss-Seidel Iteration
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                diag=eta*varkappa+alpha_diss*pow(dx,2)+6*young_modulus(rho_old[i][j]);

                u1ip1j=u1_iter_old[per_index(i+1)][j]+(per_index(i+1)+j*M<i+j*M)*(u1[per_index(i+1)][j]-u1_iter_old[per_index(i+1)][j]);
                u1im1j=u1_iter_old[per_index(i-1)][j]+(per_index(i-1)+j*M<i+j*M)*(u1[per_index(i-1)][j]-u1_iter_old[per_index(i-1)][j]);
                u1ijp1=u1_iter_old[i][per_index(j+1)]+(i+per_index(j+1)*M<i+j*M)*(u1[i][per_index(j+1)]-u1_iter_old[i][per_index(j+1)]);
                u1ijm1=u1_iter_old[i][per_index(j-1)]+(i+per_index(j-1)*M<i+j*M)*(u1[i][per_index(j-1)]-u1_iter_old[i][per_index(j-1)]);

                u1ip1jp1=u1_iter_old[per_index(i+1)][per_index(j+1)]+(per_index(i+1)+per_index(j+1)*M<i+j*M)*(u1[per_index(i+1)][per_index(j+1)]-u1_iter_old[per_index(i+1)][per_index(j+1)]);
                u1im1jm1=u1_iter_old[per_index(i-1)][per_index(j-1)]+(per_index(i-1)+per_index(j-1)*M<i+j*M)*(u1[per_index(i-1)][per_index(j-1)]-u1_iter_old[per_index(i-1)][per_index(j-1)]);
                u1ip1jm1=u1_iter_old[per_index(i+1)][per_index(j-1)]+(per_index(i+1)+per_index(j-1)*M<i+j*M)*(u1[per_index(i+1)][per_index(j-1)]-u1_iter_old[per_index(i+1)][per_index(j-1)]);
                u1im1jp1=u1_iter_old[per_index(i-1)][per_index(j+1)]+(per_index(i-1)+per_index(j+1)*M<i+j*M)*(u1[per_index(i-1)][per_index(j+1)]-u1_iter_old[per_index(i-1)][per_index(j+1)]);

                u2ip1j=u2_iter_old[per_index(i+1)][j]+(per_index(i+1)+j*M<i+j*M)*(u2[per_index(i+1)][j]-u2_iter_old[per_index(i+1)][j]);
                u2im1j=u2_iter_old[per_index(i-1)][j]+(per_index(i-1)+j*M<i+j*M)*(u2[per_index(i-1)][j]-u2_iter_old[per_index(i-1)][j]);
                u2ijp1=u2_iter_old[i][per_index(j+1)]+(i+per_index(j+1)*M<i+j*M)*(u2[i][per_index(j+1)]-u2_iter_old[i][per_index(j+1)]);
                u2ijm1=u2_iter_old[i][per_index(j-1)]+(i+per_index(j-1)*M<i+j*M)*(u2[i][per_index(j-1)]-u2_iter_old[i][per_index(j-1)]);

                u2ip1jp1=u2_iter_old[per_index(i+1)][per_index(j+1)]+(per_index(i+1)+per_index(j+1)*M<i+j*M)*(u2[per_index(i+1)][per_index(j+1)]-u2_iter_old[per_index(i+1)][per_index(j+1)]);
                u2im1jm1=u2_iter_old[per_index(i-1)][per_index(j-1)]+(per_index(i-1)+per_index(j-1)*M<i+j*M)*(u2[per_index(i-1)][per_index(j-1)]-u2_iter_old[per_index(i-1)][per_index(j-1)]);
                u2ip1jm1=u2_iter_old[per_index(i+1)][per_index(j-1)]+(per_index(i+1)+per_index(j-1)*M<i+j*M)*(u2[per_index(i+1)][per_index(j-1)]-u2_iter_old[per_index(i+1)][per_index(j-1)]);
                u2im1jp1=u2_iter_old[per_index(i-1)][per_index(j+1)]+(per_index(i-1)+per_index(j+1)*M<i+j*M)*(u2[per_index(i-1)][per_index(j+1)]-u2_iter_old[per_index(i-1)][per_index(j+1)]);

                u1[i][j]=b1_seidel[j*M+i]+young_modulus(rho_old[i][j])*((u1ijm1+u1ijp1)+2.0*(u1im1j+u1ip1j));
                u1[i][j]=u1[i][j]+young_modulus(rho_old[i][j])/4.0*(u2ip1jp1-u2ip1jm1-u2im1jp1+u2im1jm1);
                u1[i][j]=u1[i][j]+E1/4.0*(2.0*(u1ip1j-u1im1j)*(rho_old[per_index(i+1)][j]-rho_old[per_index(i-1)][j])+(u1ijp1-u1ijm1+u2ip1j-u2im1j)*(rho_old[i][per_index(j+1)]-rho_old[i][per_index(j-1)]));
                u1[i][j]=u1[i][j]/diag;

                u2[i][j]=b2_seidel[j*M+i]+young_modulus(rho_old[i][j])*(2.0*(u2ijm1+u2ijp1)+(u2im1j+u2ip1j));
                u2[i][j]=u2[i][j]+young_modulus(rho_old[i][j])/4.0*(u1ip1jp1-u1ip1jm1-u1im1jp1+u1im1jm1);
                u2[i][j]=u2[i][j]+E1/4.0*(2.0*(u2ijp1-u2ijm1)*(rho_old[i][per_index(j+1)]-rho_old[i][per_index(j-1)])+(u1ijp1-u1ijm1+u2ip1j-u2im1j)*(rho_old[per_index(i+1)][j]-rho_old[per_index(i-1)][j]));
                u2[i][j]=u2[i][j]/diag;


            }
        }



    }
    // end of Seidel for displacement


    //Seidel for chemical
    // The system will be rewritten as Ax=b
    // Compute b
    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            b_chemical[j*M+i]=varkappa*chemical[i][j]+C_0*load_old[i][j]*pow(dx,2);
        }
    }

    //iteration cycle
    seidel_c_iteration=0;
    error_c_norm=0.5;

    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            chemical_iter_old[i][j]=1.0;
        }
    }

     while ((seidel_c_iteration<1000)&&(error_c_norm>0.0001))
    {
        //count number of iterations
        seidel_c_iteration++;

        //count error_norm
        error_c_norm=0.0;

        max_c_kmkp1=0.0;
        max_c=0.0;
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                helper_seidel=fabs(chemical[i][j]-chemical_iter_old[i][j]);
                helper_seidel_2=fabs(chemical[i][j]);
                if (helper_seidel>max_c_kmkp1) max_c_kmkp1=helper_seidel;
                if (helper_seidel_2>max_c) max_c=helper_seidel_2;
            }
        }
        error_c_norm=helper_seidel/helper_seidel_2;
       // error_norm=1.0*(seidel_iteration<20);
        //old -> to -> new
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                chemical_iter_old[i][j]=chemical[i][j];
            }
        }

        //Gauss-Seidel Iteration
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                diag_c=varkappa+4*D_C+(gamma_1*rho_old[i][j]+beta_0)*pow(dx,2);

                cip1j=chemical_iter_old[per_index(i+1)][j]+(per_index(i+1)+j*M<i+j*M)*(chemical[per_index(i+1)][j]-chemical_iter_old[per_index(i+1)][j]);
                cim1j=chemical_iter_old[per_index(i-1)][j]+(per_index(i-1)+j*M<i+j*M)*(chemical[per_index(i-1)][j]-chemical_iter_old[per_index(i-1)][j]);
                cijp1=chemical_iter_old[i][per_index(j+1)]+(i+per_index(j+1)*M<i+j*M)*(chemical[i][per_index(j+1)]-chemical_iter_old[i][per_index(j+1)]);
                cijm1=chemical_iter_old[i][per_index(j-1)]+(i+per_index(j-1)*M<i+j*M)*(chemical[i][per_index(j-1)]-chemical_iter_old[i][per_index(j-1)]);

                cip1jp1=chemical_iter_old[per_index(i+1)][per_index(j+1)]+(per_index(i+1)+per_index(j+1)*M<i+j*M)*(chemical[per_index(i+1)][per_index(j+1)]-chemical_iter_old[per_index(i+1)][per_index(j+1)]);
                cim1jm1=chemical_iter_old[per_index(i-1)][per_index(j-1)]+(per_index(i-1)+per_index(j-1)*M<i+j*M)*(chemical[per_index(i-1)][per_index(j-1)]-chemical_iter_old[per_index(i-1)][per_index(j-1)]);
                cip1jm1=chemical_iter_old[per_index(i+1)][per_index(j-1)]+(per_index(i+1)+per_index(j-1)*M<i+j*M)*(chemical[per_index(i+1)][per_index(j-1)]-chemical_iter_old[per_index(i+1)][per_index(j-1)]);
                cim1jp1=chemical_iter_old[per_index(i-1)][per_index(j+1)]+(per_index(i-1)+per_index(j+1)*M<i+j*M)*(chemical[per_index(i-1)][per_index(j+1)]-chemical_iter_old[per_index(i-1)][per_index(j+1)]);

                chemical[i][j]=b_chemical[j*M+i]+D_C*(cip1j+cim1j+cijp1+cijm1);
                chemical[i][j]=chemical[i][j]/diag_c;
            }
        }
    }
    // end of Seidel for chemical


    for (i=0;i<M;i++)
    {
        for (j=0;j<M;j++)
        {
            magnitude_of_displacement[i][j]=sqrt(u1[i][j]*u1[i][j]+u2[i][j]*u2[i][j]);
            // update density
            rho[i][j]=rho_old[i][j]+dt*(-ogamma_2*chemical[i][j]*rho_old[i][j]);
        }
    }





// Updating x,y,phi
    for (n=0;n<Nc;n++)
    {
        //interaction term
        x_interaction_nk=0.0;
        y_interaction_nk=0.0;
        for (k=0;k<Nc;k++)
        {
            if (n!=k)
            {
                r_x_per=per4(x_old[k]-x_old[n]);
                r_y_per=per4(y_old[k]-y_old[n]);
                r_distance=sqrt(pow(r_x_per,2)+pow(r_y_per,2));
                force_adhesion_x=0.0;
                force_adhesion_y=0.0;
                if ((cell_label[n]==1)&&(cell_label[k]==1)) {helper_adhesion=normal_adhesion;}
                if ((cell_label[n]==2)&&(cell_label[k]==1)) {helper_adhesion=intermediate_adhesion;}
                if ((cell_label[n]==1)&&(cell_label[k]==2)) {helper_adhesion=intermediate_adhesion;}
                if ((cell_label[n]==2)&&(cell_label[k]==2)) {helper_adhesion=weak_adhesion;}

                if ((r_distance>1.0)&&(r_distance<distance_adhesion))
                {
                    force_adhesion_x=force_adhesion_x+helper_adhesion*r_x_per/r_distance;
                    force_adhesion_y=force_adhesion_y+helper_adhesion*r_y_per/r_distance;

                }
                x_interaction_nk=x_interaction_nk + cell_interaction(x_old[n],y_old[n],x_old[k],y_old[k],1)+force_adhesion_x;
                y_interaction_nk=y_interaction_nk + cell_interaction(x_old[n],y_old[n],x_old[k],y_old[k],2)+force_adhesion_y;
            }
        }


        //computation of equating displacement torque
        x_shoulder=x_old[n]-ell_shoulder*sin(phi_old[n]);
        y_shoulder=y_old[n]+ell_shoulder*cos(phi_old[n]);
        x_shoulder=per3(x_shoulder);
        y_shoulder=per3(y_shoulder);

        i = (int)((x_shoulder)/dx);
        i = min_int(i,M-1);
        i = max_int(i,0);
        j = (int)((y_shoulder)/dx);
        j = min_int(j,M-1);
        j = max_int(j,0);

        f_11=magnitude_of_displacement[i][j];
        f_12=magnitude_of_displacement[per_index(i+1)][j];
        f_21=magnitude_of_displacement[i][per_index(j+1)];
        f_22=magnitude_of_displacement[per_index(i+1)][per_index(j+1)];

        u_left[n]=interpol(x_shoulder,y_shoulder,i*dx,j*dx,f_11,f_12,f_21,f_22);

        x_shoulder=x_old[n]+ell_shoulder*sin(phi_old[n]);
        y_shoulder=y_old[n]-ell_shoulder*cos(phi_old[n]);

        x_shoulder=per3(x_shoulder);
        y_shoulder=per3(y_shoulder);


        i = (int)((x_shoulder)/dx);
        i = min_int(i,M-1);
        i = max_int(i,0);
        j = (int)((y_shoulder)/dx);
        j = min_int(j,M-1);
        j = max_int(j,0);

        f_11=magnitude_of_displacement[i][j];
        f_12=magnitude_of_displacement[per_index(i+1)][j];
        f_21=magnitude_of_displacement[i][per_index(j+1)];
        f_22=magnitude_of_displacement[per_index(i+1)][per_index(j+1)];

        u_right[n]=interpol(x_shoulder,y_shoulder,i*dx,j*dx,f_11,f_12,f_21,f_22);

        //moving on surface
        i = (int)((x_old[n])/dx);
        i = min_int(i,M-1);
        i = max_int(i,0);
        j = (int)((y_old[n])/dx);
        j = min_int(j,M-1);
        j = max_int(j,0);


        //finding spatial and angular velocities
        v1[n]=v_prop*(0.8*(cell_label[n]==1)+0.2)*rho_old[i][j]*cos(phi_old[n])+k_int*x_interaction_nk;//
        x[n]=x_old[n]+dt*v1[n];
        x[n]=per3(x[n]);

        v2[n]=v_prop*(0.8*(cell_label[n]==1)+0.2)*rho_old[i][j]*sin(phi_old[n])+k_int*y_interaction_nk;//*(0.8*(cell_label[n]==1)+0.2)
        y[n]=y_old[n]+dt*v2[n];//
        y[n]=per3(y[n]);

        ran=((double) rand() / (RAND_MAX));
        ran=2.0*ran-1.0;
        phi[n]=phi_old[n]+ran*sqrt(2*D_phi*dt);
        phi[n]=phi[n]+k_sigma*dt*(u_right[n]-u_left[n]);
    }


    //recording to files
    if (it%record_frequency==0)
    {
        //writing in console
        printf("time = %f, number of iterations = %d, i=%d, j=%d, phi/pi=%f, error=%1.10f; \n",it*dt,seidel_iteration,i,j,phi[0]/M_PI,error_norm);

        //substrate
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                fprintf(substrate,"%f ",rho[i][j]);
            }
            fprintf(substrate,"\n");
        }
        //record load
        /*
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                fprintf(loads,"%f ",load[i][j]);
            }
            fprintf(loads,"\n");
        }
        //record chemical
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                fprintf(chemical_file,"%f ",chemical[i][j]);
            }
            fprintf(chemical_file,"\n");
        }*/

        //record magnitude of displacement
        for (i=0;i<M;i++)
        {
            for (j=0;j<M;j++)
            {
                fprintf(displacement_magnitude,"%f ",magnitude_of_displacement[i][j]);
            }
            fprintf(displacement_magnitude,"\n");
        }
        //cell trajectories
        for (n=0;n<Nc;n++)
        {
            fprintf(cells,"%f %f %f %f %f %f %f %d\n",x[n],y[n],v1[n],v2[n],phi[n],u_left[n],u_right[n],cell_label[n]);
        }

    }
}


//Close all edited files
fclose(substrate);
fclose(cells);
fclose(loads);
fclose(chemical_file);
fclose(displacement_magnitude);

}

return 0;
}
