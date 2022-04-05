#include <iostream>
#include <math.h>
#include <fstream>

using namespace std;

int M1 = 10;
int DT = 2;
int ALPHA = 0.0000001;
int DX2 = 0.000016;
int NTST = 50;

void init_phi(double* phi)
{
    for(int i=1;i<=M1;i++)
        phi[i]=0.0;
}

void bcond(double* phi)
{
    phi[0]=-100.0;
    //phi[M1+1]=200.0;
    phi[M1+1]=phi[M1];
}

void calc_source(double* phi,double* source)
{
    for(int i=1;i<=M1;i++)
    {
        source[i]=1.0/DT*phi[i]+0.5*ALPHA*(phi[i+1]-2.*phi[i]+phi[i-1])/DX2;
    }
    source[1]+=0.5*ALPHA*phi[0]/DX2;
}

void calc_coef(double* coef, double* source, const int& nTime)
{
    double coefRef[M1+1][M1+1];
    for(int i=1;i<=M1;i++)
        for(int j=1;j<=M1;j++)
            coefRef[i][j]=0.0;
    for(int i=1;i<=M1;i++)
    {
        coefRef[i][i]=1.0/DT+0.5*ALPHA*2.0/DX2;
        if(i!=1) coefRef[i-1][i]=-0.5*ALPHA/DX2;
        if(i!=M1) coefRef[i+1][i]=-0.5*ALPHA/DX2;
    }
    coefRef[M1][M1]=1.0/DT+0.5*ALPHA*1.0/DX2;
    for(int i=1;i<=M1;i++)
        for(int j=1;j<=M1;j++)
            coef[i*M1+j]=coefRef[i][j];
}

void solve(double* coef, double* phi, double* source, const int& nTime)
{
    int nIter=1;
    double coefRef[M1+1][M1+1], oldPhi[M1+2], residual;

    for(int i=1;i<=M1;i++)
        for(int j=1;j<=M1;j++)
            coefRef[i][j] = coef[i*M1+j];
    do
    {
        for(int i=1;i<M1;i++)
            oldPhi[i]=phi[i];
        for(int i=1;i<M1;i++)
        {
            double sumLHS = 0.0;
            for(int iLHS=1;iLHS<=M1;iLHS++)
                sumLHS+=coefRef[iLHS][i]*oldPhi[iLHS];
            sumLHS-=coefRef[i][i]*oldPhi[i];
            phi[i]=(source[i]-sumLHS)/(coefRef[i][i]);
        }
        residual=0.0;
        for(int i=1;i<=M1;i++)
            residual+=fabs(phi[i]-oldPhi[i])/double(M1);
        cout << nTime << " : RES = " << residual << ", ITER= " << nIter << endl;
        nIter++;
    } while (residual>0.0001);
}

void write_phi(double* phi)
{
    ofstream os("1D.dat");
    for(int i=1;i<=M1;i++)
        os << phi[i] << endl;
    os.close();
}

int main()
{
    double phi[M1+2], source[M1+1], coef[M1+1][M1+1];
    init_phi(phi);
    for(int nTime=0;nTime<3;nTime++)
    {
        bcond(phi);
        cout << phi[0] << " : " << phi[M1] << endl;
        calc_source(phi,source);
        // calc_coef((double*)coef,source,nTime);
        // solve((double*)coef,phi,source,nTime);
    }
    // write_phi(phi);
    for(int i=1;i<=M1;i++)
        cout << source[i]<<" : ";
    //1.0/DT*phi[i]+0.5*ALPHA*(phi[i+1]-2.*phi[i]+phi[i-1])/DX2
    for(int i=1;i<=M1;i++)
        for(int j=1;j<=M1;j++)
            cout << coef[i][j]<<" : ";
    //cout << coef << " : " << phi[1] << " : " << phi[2] << endl;
    return 0;                                  
}