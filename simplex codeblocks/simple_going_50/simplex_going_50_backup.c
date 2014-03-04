//#include<iostream>
#include <cstdlib>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/core>
#include<cfloat>
using namespace std;
using namespace Eigen;

int n,m;
MatrixXd degenerate(MatrixXd b, long double episol)
{
    int i;
    for(i=0;i<m;i++)
        {
            b(i,0)=b(i,0)+pow(episol,i+1);
        }
    return b;
}
MatrixXd convert_b(MatrixXd b, MatrixXd dummy_b)
{
    int i,j;
    for(i=0;i<m-n;i++)
    {
        b(i,0)=dummy_b(i,0);
    }
    for(;i<m;i++)
        b(i,0)=0;
    return b;
}

MatrixXd convert_c(MatrixXd c, MatrixXd dummy_c)
{
    int i,j;
    for(i=0,j=0;i<n/2;i++)
    {
        c(j,0)=dummy_c(i,0);
        c(j+1,0)=-dummy_c(i,0);
        j+=2;
    }
    return c;
}
MatrixXd convert_A(MatrixXd A, MatrixXd dummy_a)
{
    int i, j,k=0,l,h;
    for(i=0;i<m-n;i++)
        for(j=0,l=0,h=0;j<n/2;j++)
        {
            A(i,h)=dummy_a(i,j);
            A(i,h+1)=-dummy_a(i,j);
            h+=2;
        }
    for(;i<m;i++,k++){
        for(j=0;j<n;j++)
        {
            A(i,k)=-1;
            if(j!=k)
                A(i,j)=0;
        }
    }
return A;
}

bool System_unbounded( const MatrixXd& A2,const MatrixXd& B1)
{
    for(int l=0;l<A2.rows();l++)
        if((A2.row(l)*B1)(0)>0)
          return false;
    return true;
}

MatrixXd get_solution(MatrixXd Ax,MatrixXd b,MatrixXd ct,MatrixXd A,MatrixXd x, MatrixXd Bt,MatrixXd B1){
    int count;
    int dege=0;
    int first_run=1;
    long double episol = 0.1;//pow(16,1.0/m);
    MatrixXd b_value=b;
   // bool non_degenerate = false;
while(1){
    cout<<"===================================="<<endl;
    cout<<"x="<<x<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;

    Ax=A*x;
    cout<<"====================================="<<endl;
    cout<<"Ax="<<Ax<<endl;
    cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;
    count=0;
    float g,h;
    for(int i=0;i<m;i++)
    {
        g=Ax(i,0);
        h=b(i,0);
        if(abs(g-h)<0.00001)
        {

            count++;
            cout<<"====================================="<<endl;
            cout<<"Ax="<<Ax(i,0)<<endl;
            cout<<"b="<<b(i,0)<<endl;
            cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;

        }
        else
        {

            cout<<"====================================="<<endl;
            cout<<"Ax="<<Ax(i,0)<<endl;
            cout<<"b="<<b(i,0)<<endl;
            cout<<"+++++++++++++++++++++++++++++++++++++"<<endl;

        }
    }
    cout<<"-----------------"<<endl;
    cout<<count<<endl;
    cout<<"-------------------"<<endl;
    MatrixXd A1(count,n);
    MatrixXd A2((m-count), n);
    int count1=0;
    int count2=0;
    int dummy[A2.rows()];
    for(int i=0;i<m;i++)    //creating my A1 and A2 matrices
    {
        g=Ax(i,0);
        h=b(i,0);
        if(abs(g-h)<0.00001)
        {
            for(int j=0;j<n;j++){
                A1(count1,j)=A(i,j);
            }
            count1++;
        }
        else
        {
            for(int j=0;j<n;j++)
            {
                A2(count2,j)=A(i,j);
            }
            dummy[count2]=i;
            count2++;
        }
    }
    cout<<" Hetre!!! The Matrices are ::"<<A1<<" and ::"<<endl<<A2<<endl;

cout<<"========================================"<<endl;
if(count1>n){
    dege++;
    if(dege!=1)
        episol=episol/2;
    b=degenerate(b_value,episol);
    cout<<"b=\n"<<b<<endl;
    continue;
}
else if(count1<n){

    episol=episol/2;
    b=degenerate(b_value,episol);
    cout<<"b=\n"<<b<<endl;
    continue;
}
cout<<"============================================="<<endl;
MatrixXd B;
B= -1*(A1.inverse());
cout<<"================================"<<endl;
cout<<B<<endl;
cout<<"================================"<<endl;
float k=0;
int rows_a2[A2.rows()];
float av[A2.rows()];
for(int i=0;i<A2.rows();i++){
    rows_a2[i]=0;
    av[i]=0;
}
int z=0;
int stat=0;
for(int i=0;i<B.cols();i++)
{
    k =0;
    Bt=B.col(i);
    cout<<"==============================="<<endl;
    cout<<"Bt="<<Bt<<endl;
    cout<<"+++++++++++++++++++++++++++++++"<<endl;
    if((ct*Bt)(0)>0)
    {
        stat++;
        B1=Bt;
        for(int j=0;j<A2.rows();j++)
        {
            float sum=0;
            for(int k=0;k<n;k++)
            {
                sum+=A2(j,k)*B1(k,0);
            }

            if(sum>=0)
            {
                cout<<"sum="<<sum<<endl;
                av[z]=sum;
                rows_a2[z]=j;
                z++;              //z defines how many rows of A2 are there such that A2*bt>0
            }
        }
        break;
    }
}
if(stat==0)
    break;
float t;
int u;
float ax=0,mini;
int one=1;
if(System_unbounded(A2,B1))
{
    cout<<"System is unbounded"<<endl;
    exit(1);
}
for(int a=0;a<z;a++)
{
    ax = 0;

    ax = (A2.row((rows_a2[a]))*x)(0);
    for(int l=0;l<A.rows();l++)
    {
        if(A2.row((rows_a2[a]))== A.row(l))
            u=l;
    }

    cout<<"ax="<<ax<<endl;
    t= ((b(u,0)-ax)/av[a]);
    cout<<"t="<<t<<endl;
    if(one==1)
    {
        mini=t;
        one++;
    }
    if(t<mini)
        mini=t;
    cout<<"mini="<<mini<<endl;
}   x=x+mini*B1;
cout<<x<<endl;
cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<mini<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    //break;
}
return x;


}

MatrixXd get_vertex(const MatrixXd& a, const MatrixXd& b, MatrixXd& x)
{
   int count;
   MatrixXd vert(n,1);
   MatrixXd x_dummy(n,1);
   x_dummy=x;
   bool got_vertex=false;
   MatrixXd ax(m,1);
   while(!got_vertex){
    ax=a*x_dummy;
    cout<<"\n Ax= "<<ax<<endl;
    count=0;
    float g,h;
    for(int i=0;i<m;i++)
    {
        g=ax(i,0);
        h=b(i,0);
        if(abs(g-h)<0.00001)
        {
            count++;
        }
        else
        {}    //nothing to do
    }
    cout<<"count=\n"<<count<<endl;
    MatrixXd A1(count,n);
    MatrixXd A2((m-count), n);
    int count1=0;
    int count2=0;
    int dummy[A2.rows()];
    for(int i=0;i<m;i++)    //creating my A1 and A2 matrices
    {
        g=ax(i,0); h=b(i,0);
        if(abs(g-h)<0.00001)
        {
            for(int j=0;j<n;j++){
                A1(count1,j)=a(i,j);
            }
            count1++;
        }
        else
        {
            for(int j=0;j<n;j++)
            {
                A2(count2,j)=a(i,j);
            }
            dummy[count2]=i;
            count2++;
        }
    }
    cout<<" Matrices are A1::\n"<<A1<<" \n A2::"<<endl<<A2<<endl;
    MatrixXd B;
    int z,run=1;
 if(A1.rows()!=0){
    FullPivLU<MatrixXd> lu_decomp(A1);
    if(lu_decomp.rank()==n)
    {
        cout<<"INTOTHEBOKERS\n"<<endl;
        cout<<"A1=\n"<<A1<<endl;
        got_vertex=true;
        continue;
    }
    vert=(lu_decomp.kernel()).col(0);
 }
 else
 {
    for(int j=0;j<n;j++)
        vert(j,0)=0.1;
 }
bool flag=false;  //edited if beta negative
    if(System_unbounded(A2,vert))
        flag=true;//edited if
 if(vert.maxCoeff()==0 && vert.minCoeff()==0)
 {
     cout<<"vertex is zero"<<endl;
    got_vertex=true;
    continue;
 }
 cout<<"vertex=\n"<<vert<<endl;
 float beta,test;
 for(int i=0;i<A2.rows();i++)
 {
    for(int l=0;l<a.rows();l++)
    {
        if(A2.row(i)== a.row(l))
            z=l;
    }
    cout<<"b="<<b(z,0)<<"A2.row="<<A2.row(i)*x_dummy<<"A2.rowsssssssssss="<<A2.row(i)*vert<<endl;
     test= (b(z,0)-(A2.row(i)*x_dummy)(0))/(A2.row(i)*vert)(0);
    cout<<"test="<<test<<endl;
    if(flag==false){                   //edited if beta neg.
    if(run==1 && test>0){
        beta=test;
        run++;
    }
    if(beta>test && test>0)
        beta=test;
    }
    else     {            //Added to take into account if all beta's are negative and then move in -beta direction
        if(run==1){
        beta=test;
        run++;
    }
    if(beta<test)
        beta=test;
    }

 }
 cout<<"beta="<<beta;
 x_dummy=x_dummy+beta*vert;
 cout<<"x_dummy=\n"<<x_dummy<<endl;
   }

   cout<<"\nFinal Intial vertex=\n"<<x_dummy<<endl;
   cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
   return x_dummy;
}




int main(int argc, char** argv) {
    cout <<"Enter the no. of Rows and Columns of A:";
    cin>>m>>n;
    MatrixXd A(m,n),b(m,1),c(n,1),x(n,1),Ax(m,1),B1(n,1),Bt(n,1),ct(1,n),vertex(n,1),sol(m,n);
    bool rank_prb=false; //mxn matrix
    cout <<"Enter the 'A' matrix row-wise:"<<endl;
    for(int i=0;i<m;i++)
        for(int j=0;j<n;j++)
            cin>>A(i,j);

    cout <<"Enter the 'b' matrix:"<<endl;
    for(int j=0;j<m;j++)
        cin>>b(j,0);

    cout <<"Enter the 'c' cost matrix:"<<endl;
    for(int k=0;k<n;k++)
        cin>>c(k,0);


    cout <<"Enter the initial feasible solution 'x':"<<endl;

    for(int l=0;l<n;l++)
        cin>>x(l,0);
    FullPivLU<MatrixXd> lu_decomp(A);
    if(lu_decomp.rank()!=n){
        MatrixXd dummy_a, dummy_b, dummy_c,dummy_x;
        dummy_a=A;dummy_b=b;dummy_c=c;dummy_x=x;
        n=2*n;m=m+n;
        A.resize(m,n);b.resize(m,1);c.resize(n,1);x.resize(n,1);
        Ax.resize(m,1);ct.resize(1,n);Bt.resize(n,1);B1.resize(n,1);
        A=convert_A(A,dummy_a);
        cout<<"Done A\n";
        b=convert_b(b,dummy_b);
        cout<<"Done b\n";
        c=convert_c(c,dummy_c);
        x=convert_c(x,dummy_x);

        rank_prb=true;
    }
    cout<<"A=\n"<<A<<endl;
    cout<<"b=\n"<<b<<endl;
    cout<<"c=\n"<<c<<endl;
    ct=c.transpose();
    x=get_vertex(A,b,x);
    x=get_solution(Ax,b,ct,A,x,Bt,B1);
    if(!rank_prb)
        sol=x;
    else{
        int i,k;
        {
            sol.col(0)=x.col(0);
        }

    }
cout<<"Optimum solution is "<<endl;
cout <<sol<<endl;
cout<<"Cost is:"<<ct*sol<<endl;

    return 0;
}

