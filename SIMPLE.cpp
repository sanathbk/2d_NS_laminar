#include<iostream>
#include<fstream>
#include<cmath>
#include <sys/types.h>
#include <sys/stat.h>

using namespace std;

const double Re=50.0 ; //Reynolds Number
const double length1= 12.0;//Length
const double length2= 1.0;//Height
const double dt= 1;// TimeStep
const int N= 201;
const int M= 51;
const double rho= 1.0;
double tolerance=pow(10,-6);
const double NT=300000;//Max number of iterations allowed

// defining variables
double u_infinity= 1.0;

double p[N+1][M+1],p_error[N+1][M+1]={0.0};
double u[N+1][M+1], un[N+1][M+1],u_n_1[N+1][M+1]={0.0};
double v[N+1][M+1], vn[N][M+1],v_n_1[N+1][M+1]={0.0};
double l[N+1],b[N+1],c[N+1],d[N+1],x[N+1], xm[M+1];//TDMA Coefficients

//Momentum Equation Coefficient Functions
double a(int ix, int iy,int z);
double a_w(int ix, int iy,int z);
double a_e(int ix, int iy,int z);
double a_s(int ix, int iy,int z);
double a_n(int ix, int iy,int z);

//Pressure Correction Equation Coefficients
double p_e(int ix, int iy);
double p_s(int ix, int iy);
double p_n(int ix, int iy);
double p_w(int ix, int iy);
double dx= length1/(N-1);
double dy= length2/(M-1);

//double viscosity= (rho*u_infinity*length2)/rey ;
//double D1= viscosity/dx;
//double D2= viscosity/dy;

double R1=dy/dx;
double R2=dx/dy;

//Momentum-Equation Solver
void gs_u(void);
void gs_v(void);

//Pressure Correction Solver
void gs_p(void);

void TDMAm(void);//TDMA X-Sweep_Momentum Eqns
void TDMAn(void);//TDMA Y-Sweep_Momentum Eqns
void TDMApm(void);


int main()
	{
		if (mkdir("./TimeData", ACCESSPERMS)!=0);
			cout<<"Couldn't Create Directory"<<endl;
        double t = 0;
        int j = 0;
        double maxu = 0, erru = 0;
        double max = 1;
        ofstream fout20;
        fout20.open("Convergence Plot.dat");

        while(max > 1e-6)//1e-6 is convergence criteria for steady state

        {        	
            maxu = 0;
            j++;
            t =j*dt;
            cout<<"Time = "<<t<<"s"<<endl<<endl;
            char time_data[50];
            ofstream fout15;
            fout15.open("tolerance.dat");

        	// Updating previous timestep values-u velocity
        	for(int ix=0;ix<N;ix++)
        	{
            	for(int iy=0;iy<=M;iy++)
            	{
                	u_n_1[ix][iy]=u[ix][iy];//u_n_1 is velocity of previous timestep
            	}
       		}

       		// Updating previous timestep values-v velocity
       		for(int ix=0;ix<=N;ix++)
       		{
        		for(int iy=0;iy<M;iy++)
        		{
            		v_n_1[ix][iy]=v[ix][iy];
        		}
       		}

////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
 //                       SIMPLE algorithm
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

	 		for(int it=NT;it>0;it--)
    		{

        		cout<<"Iteration number: "<<NT- it<<endl;

          		//Boundary Conditions-Velocity

    			for (int i=1;i<M;i++)
    			{
        			u[0][i]=1.0;//Left_Velocity Inlet
        			u[N-1][i]=u[N-2][i];//Right
        		}
        		for (int i=1;i<N;i++)
        		{		
        			v[i][0]=0.0;//Bottom
        			v[i][M-1]=0.0;//Top
        		}
    			for (int i=0;i<=M-1;i++)
    			{
        			v[0][i]=-v[1][i];//Left
        			v[N][i]=-v[N-1][i];//Right
        		}
    			for (int i=0;i<=N-1;i++)
    			{
        			u[i][0]=-u[i][1];//Bottom
        			u[i][M]=-u[i][M-1];//Top
        		}

        		//Boundary Conditions-Pressure
         		for(int ix=0;ix<=M;ix++)
         		{
            		p[N][ix]=p[N-1][ix];//Right
        		}
        		for(int ix=0;ix<=N;ix++)
        		{
            		p[ix][M]=p[ix][M-1];//Top
        		}
         		for(int ix=0;ix<=M;ix++)
         		{
            		p[0][ix]=p[1][ix];//Left
        		}
        		for(int ix=0;ix<=N;ix++)
        		{
            		p[ix][0]=p[ix][1];//Bottom
        		}
				//p[0][0]=0.0;

				//copying values of U-velocity of previous iteration
        		for(int ix=0;ix<N;ix++)
        		{
            		for(int iy=0;iy<=M;iy++)
            		{
                		un[ix][iy]=u[ix][iy];
            		}
       			}
				//copying values of V-velocity of previous iteration

       			for(int ix=0;ix<=N;ix++)
       			{
        			for(int iy=0;iy<M;iy++)
        			{
            			vn[ix][iy]=v[ix][iy];
        			}
       			}


        		// solving the system of linear equations for u and v momentum equations and pressure correction equation
        		gs_u();
        		gs_v();

        		for(int ix=0;ix<N;ix++)
        		{
            		for(int iy=0;iy<M;iy++)
            		{
                		p_error[ix][iy]=0.0;
            		}
       			}

       			//Solving Pressure Correction Equation
        		gs_p();

          		//updating pressure with under-relaxation factor

            	for(int ix=1;ix<N;ix++)
            	{
                	for(int iy=1;iy<M;iy++)
                	{
                    	p[ix][iy]= p[ix][iy]+ 0.7 * p_error[ix][iy];// here '0.7' is the under relaxation  factor
                   	}
            	}

            	//update velocity values using correction equations

            	//updating u velocity

            	for(int ix=N-2;ix>0;ix--)
            	{
                	for(int iy=M-1;iy>0;iy--)
                	{
                        u[ix][iy]= u[ix][iy]+ (dy*1.0/a(ix,iy,1))*(p_error[ix][iy]-p_error[ix+1][iy]);//updating velocity values without using under relaxation factor
                    }
                }
                //updating v velocity
            	for(int ix=N-1;ix>0;ix--)
            	{
                	for(int iy=M-2;iy>0;iy--)
                	{
                		v[ix][iy]= v[ix][iy]+ (dx*1.0/a(ix,iy,2))*(p_error[ix][iy]-p_error[ix][iy+1]);//updating velocity values without using under relaxation factor
                	}
            	}
            	//Printing values to check if code is working
            	cout<<u[N/2][M/2]<<endl;
            	cout<<u[N][M/2]<<endl;


            	//checking convergence

            	double maximum1=0.0;
            	double maximum2=0.0;
            	double maximum=0.0;

            	for(int ix=0;ix<N-1;ix++)
            	{
                	for(int iy=0;iy<M;iy++)
                	{
    	                if(abs(u[ix][iy]-un[ix][iy])>maximum1)
    	                {
                        	maximum1= abs(u[ix][iy]-un[ix][iy]);//u is velocity in current iteration; un is velocity of previous iteration
                       	}
         	       	}
            	}

            	for(int ix=0;ix<N;ix++)
            	{
                	for(int iy=0;iy<M-1;iy++)
                	{
    		           	if(abs(v[ix][iy]-vn[ix][iy])>maximum1)
    		           	{
                        	maximum1= abs(v[ix][iy]-vn[ix][iy]);
                        }
                	}
            	}

            	if(maximum1<maximum2)
            	{
                	maximum = maximum2;
            	}
            	else if(maximum2<maximum1)
            	{
                	maximum = maximum1;
            	}
            	if(maximum<tolerance)
                break;
    			cout<<"Value of error:"<<maximum<<endl;//Printing value of error after each iteration
    			fout15<<NT-it<<" "<<maximum<<endl;
    		}

			cout<<"Writing values to file."<<endl;
           //write values to a file after one time step converges
            snprintf(time_data,sizeof time_data, "./TimeData/Data_%.2f.txt",t);
            ofstream ofp1(time_data, ios::out);
            int m;
           	for(int i=0;i<M;i++)
           	{
                m=25;
                ofp1<<0.5*(u[10][i]+u[10][i+1])/u_infinity<<" ";
                while (m<=N)
                {
                    ofp1<<(u[m][i]+u[m][i+1])*0.5/u_infinity<<" ";
                    m=m+25;
                }
                ofp1<<(dy*i)/length2<<endl;
            }

				fout15.close();
			for(int ix=0;ix<N;ix++)
			{
           		for(int iy=0;iy<M;iy++)
           		{
            		erru = abs(u[ix][iy] - u_n_1[ix][iy]);
            		if(erru > maxu)
                	maxu =erru;
           		}
       		}
       		max = maxu;
    		cout<<"Velocity Error at the end of timestep: "<<max<<endl;
    		ofstream fout25;
    		fout25.open("Velocity Contour.dat");
    		fout20<<t<<"  "<<u[N/2][M/2]<<endl;
    		for(int ix=1;ix<N+1;ix++)
			{
           		for(int iy=1;iy<M+1;iy++)
           		{
            		fout25<<t<<" "<<(ix-1)*dx<<" "<<(iy-1)*dy<<" "<<0.5*(u[ix-1][iy]+u[ix-1][iy-1])<<" "<<0.5*(v[ix-1][iy-1]+v[ix][iy-1])<<endl;
           		}
       		}
    	}
	}


// momentum equation coefficients

double a_n(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix+1][iy])-R2/Re;
     else if(z==2)
     ans= dx*0.25*(v_n_1[ix][iy]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}

double a_s(int ix, int iy,int z){

    double ans=0.0;

        if(z==1)
        ans= -dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix+1][iy-1])-R2/Re;
        else if(z==2)
        ans= -dx*0.25*(v_n_1[ix][iy-1]+v_n_1[ix][iy+1])-R2/Re;

    return ans;
}


double a_e(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix+1][iy])-R1/Re;
      else if(z==2)
      ans= dy*0.25*(u_n_1[ix][iy]+u_n_1[ix][iy+1])-R1/Re;

    return ans;
}

double a_w(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
    ans= -dy*0.25*(u_n_1[ix][iy]+u_n_1[ix-1][iy])-R1/Re;
    else if(z==2)
    ans= -dy*0.25*(u_n_1[ix-1][iy]+u_n_1[ix-1][iy+1])-R1/Re;


    return ans;
}


double a(int ix, int iy,int z){

    double ans=0.0;

    if(z==1)
     //ans= dx*dy/dt +(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
     ans= (dx*dy/dt)+(a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1))+4.0*(R1/Re+R2/Re);
      //ans= (a_w(ix,iy,1)+a_e(ix,iy,1)+a_n(ix,iy,1)+a_s(ix,iy,1));
     else if(z==2)
     //ans= dx*dy/dt +(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    ans= (dx*dy/dt)+(a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2))+4.0*(R1/Re+R2/Re);
    //ans= (a_w(ix,iy,2)+a_e(ix,iy,2)+a_n(ix,iy,2)+a_s(ix,iy,2));
    return ans;
}

// pressure correction equation coefficients

double p_e(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix,iy,1);

    return answer;

}

double p_n(int ix, int iy){

    double answer=0.0;


    answer= -1.0*dx*dx*1.0/a(ix,iy,2);


    return answer;

}

double p_w(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dy*dy*1.0/a(ix-1,iy,1);

    return answer;

}

double p_s(int ix, int iy){

    double answer=0.0;

    answer= -1.0*dx*dx*1.0/a(ix,iy-1,2);


    return answer;

}
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void gs_u(void){
//Y-Sweep Substitution
for(int iy=1;iy<M;iy++){

    for(int ix=1;ix<N-1;ix++){

      if(ix==1)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix-1][iy]*a_w(ix,iy,1)+(dx*dy/dt)*u_n_1[ix][iy];
            l[ix]=0; c[ix]=a_e(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        if(ix==N-2)
        {
            d[ix]=(p[ix][iy]-p[ix+1][iy])*dy -u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))-u[ix+1][iy]*a_e(ix,iy,1)+(dx*dy/dt)*u_n_1[ix][iy];
            l[ix]=a_w(ix,iy,1); b[ix]=a(ix,iy,1);
            continue;
        }
        l[ix]=a_w(ix,iy,1);
        c[ix]=a_e(ix,iy,1);
        b[ix]=a(ix,iy,1);
        d[ix]=(p[ix][iy]-p[ix+1][iy])*dy-u[ix][iy+1]*(a_n(ix,iy,1))-u[ix][iy-1]*(a_s(ix,iy,1))+(dx*dy/dt)*u_n_1[ix][iy];
                  }
    TDMAn();
    for(int ix=1;ix<N-1;ix++)
        u[ix][iy]=x[ix];

}
//X-Sweep Substitution
/*for(int ix=1;ix<N-1;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy - u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy-1]*a_s(ix,iy,1);
            l[iy]=0; c[iy]=a_n(ix,iy,1); b[iy]=a(ix,iy,1);
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1))-u[ix][iy+1]*a_n(ix,iy,1);
            l[iy]=a_s(ix,iy,1); b[iy]=a(ix,iy,1);c[iy] = 0;
            continue;
        }
        l[iy]=a_s(ix,iy,1);
        //cout<<l[iy]<<endl;
        c[iy]=a_n(ix,iy,1);
        b[iy]=a(ix,iy,1);
        d[iy]=(p[ix][iy]-p[ix+1][iy])*dy- u[ix+1][iy]*(a_e(ix,iy,1))-u[ix-1][iy]*(a_w(ix,iy,1));
                  }
    TDMAm();
    for(int iy=1;iy<=M-1;iy++)
        {u[ix][iy]=xm[iy];
        //cout<<u[ix][iy]<<endl;
        }
}

*/
}

void gs_p(void){
//Y-Sweep Substitution
for(int iy=1;iy<M;iy++){

   for(int ix=1;ix<N;ix++){
     if(ix==1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix-1][iy]*p_w(ix,iy);
            l[ix]=0; c[ix]=p_e(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(ix==N-1)
        {
            d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy)-p_error[ix+1][iy]*p_e(ix,iy);
            l[ix]=p_w(ix,iy); b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[ix]=p_w(ix,iy);
        c[ix]=p_e(ix,iy);
        b[ix]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[ix]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy-p_error[ix][iy-1]*p_s(ix,iy) - p_error[ix][iy+1]*p_n(ix,iy);
        }
    TDMApm();
    for(int ix=1;ix<N;ix++)
        p_error[ix][iy]=x[ix];
}
//X-Sweep Substitution
/*
for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M;iy++){
         if(iy==1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy-1]*p_s(ix,iy);
            l[iy]=0; c[iy]=p_n(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        if(iy==M-1)
        {
            d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy)-p_error[ix][iy+1]*p_n(ix,iy);
            l[iy]=p_s(ix,iy); b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
            continue;
        }
        l[iy]=p_s(ix,iy);
        c[iy]=p_n(ix,iy);
        b[iy]=-1.0*(p_n(ix,iy)+p_w(ix,iy)+p_e(ix,iy)+p_s(ix,iy));
        d[iy]=v[ix][iy-1]*dx-v[ix][iy]*dx+u[ix-1][iy]*dy-u[ix][iy]*dy -p_error[ix-1][iy]*p_w(ix,iy) - p_error[ix+1][iy]*p_e(ix,iy);
                  }
    TDMAm();
    for(int iy=1;iy<M;iy++)
        p_error[ix][iy]=x[iy];
}*/
}


void gs_v(void){

for(int ix=1;ix<N;ix++){

    for(int iy=1;iy<M-1;iy++){

        v[ix][iy]= (dx*1.0*((p[ix][iy]-p[ix][iy+1]))-1.0*(v[ix+1][iy]*(a_e(ix,iy,2)) + v[ix][iy+1]*(a_n(ix,iy,2)) + v[ix-1][iy]*(a_w(ix,iy,2)) + v[ix][iy-1]*(a_s(ix,iy,2)))+v_n_1[ix][iy]*dx*dy/dt)/(a(ix,iy,2));

    }
}

}
//Y-Sweep Elimination
void TDMAn(void){
 int i;
 double t;
    for (i=2;i<=N-2;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-2;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
        t=x[i];
    }
}
//X-Sweep Elimination
void TDMAm(void){
 int i;
    for (i=2;i<=M-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=M-1;
    xm[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        xm[i]=(d[i]-c[i]*xm[i+1])/b[i];
    }
}

void TDMApm(void){
 int i;
    for (i=2;i<=N-1;i++)
    {
        l[i]=l[i]/(b[i-1]);
        b[i]=b[i]-l[i]*c[i-1];
        d[i]=d[i]-l[i]*d[i-1];
    }
    int n=N-1;
    x[n]=d[n]/b[n];
    for (i=n-1;i>=1;i--)
    {
        x[i]=(d[i]-c[i]*x[i+1])/b[i];
    }
}
