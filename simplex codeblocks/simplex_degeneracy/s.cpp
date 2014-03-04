#include<iostream>
#include <cstdlib>
#include<cmath>
#include <Eigen/Dense>
#include <Eigen/LU>
#include <Eigen/core>
#include<cfloat>
using namespace std;
using namespace Eigen;


MatrixXd degenerate(MatrixXd b, double episol, int m,int n)
{
    int i;
    for(i=0;i<m;i++)
        {
            b(i,0)=b(i,0)+pow(episol,i+1);
        }

    cout<<"Degenerate b="<<b<<endl;
    int dummy;
    //cin>>dummy;
    return b;
}

MatrixXd convert_x(MatrixXd x, MatrixXd dummy_x,int m ,int n )
{
    int i,j;
    for(i=0,j=0;i<n/2;i++)
    {
        if(dummy_x(i,0)>=0)
        {
            x(j,0)=dummy_x(i,0);
            x(j+1,0)=0;
            j+=2;
        }
        else
        {
            x(j,0)=0;
            x(j+1,0) = -1*dummy_x(i,0);
            j+=2;
        }
    }
    return x;
}
MatrixXd convert_b(MatrixXd b, MatrixXd dummy_b,int m,int n)
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

MatrixXd convert_c(MatrixXd c, MatrixXd dummy_c,int m,int n)
{
    int i,j;
    for(i=0,j=0;i<n/2;i++)
    {
        c(j,0)=dummy_c(i,0);
        if(dummy_c(i,0)==0)
            c(j+1,0)=dummy_c(i,0);
        else
            c(j+1,0)=-dummy_c(i,0);
        j+=2;
    }
    return c;
}
MatrixXd convert_A(MatrixXd A, MatrixXd dummy_a,int m,int n)
{
    int i, j,k=0,l,h;
    for(i=0;i<m-n;i++)
        for(j=0,l=0,h=0;j<n/2;j++)
        {
            A(i,h)=dummy_a(i,j);
            if(dummy_a(i,j)==0)
                A(i,h+1)=dummy_a(i,j);
            else
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

bool System_unbounded( const MatrixXd A2,const MatrixXd B1)
{
    for(int l=0;l<A2.rows();l++)
        if((A2.row(l)*B1)(0)>0)
          return false;
    return true;
}

MatrixXd get_vertex(const MatrixXd a, const MatrixXd b, MatrixXd x, int m ,int n)
{
   int count;
   MatrixXd vert(n,1);
   MatrixXd x_dummy(n,1);
   MatrixXd ax(m,1);
   bool got_vertex=false;
   x_dummy=x;
   while(!got_vertex){
        ax=a*x_dummy;
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
        int count1=0,count2=0;
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
                cout<<"INTOTHEBokers\n"<<endl;
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
        cout<<"direction=\n"<<vert<<endl;
        float beta,test;
        for(int i=0;i<A2.rows();i++)
        {
            for(int l=0;l<a.rows();l++)
            {
                if(A2.row(i)== a.row(l))
                    z=l;
            }
            if(!(abs((A2.row(i)*vert)(0)-0)<0.0000000001))
                test= (b(z,0)-(A2.row(i)*x_dummy)(0))/(A2.row(i)*vert)(0);
            else
                continue;
            cout<<"test="<<test<<endl;
            if(flag==false){                   //edited if beta neg.
                if(run==1 && test>0){
                    beta=test;
                    run++;
                }
                if(beta>test && test>0)
                    beta=test;
                }
            else{            //Added to take into account if all beta's are negative and then move in -beta direction
                if(run==1){
                    beta=test;
                    run++;
                }
                if(beta<test)
                    beta=test;
                }

        }
        cout<<"beta="<<beta<<endl;
        x_dummy=x_dummy+beta*vert;
        cout<<"x_dummy=\n"<<x_dummy<<endl;
    }
    cout<<"\nFinal Intial vertex=\n"<<x_dummy<<endl;
    cout<<"%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%"<<endl;
    return x_dummy;
}

MatrixXd get_solution(MatrixXd Ax,MatrixXd b,MatrixXd ct,MatrixXd A,MatrixXd x, MatrixXd Bt,MatrixXd B1,int m, int n){
    int count;
    int first_run=1;
    long double episol = 0.1;
    MatrixXd b_value=b;
    MatrixXd A_copy;
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
        static int iteration;
        iteration++;
        cout<<"iteration="<<iteration;
        cout<<"m="<<m<<endl;
        cout<<"b="<<b<<endl;
        for(int i=0;i<m;i++)
        {
            g=Ax(i,0);
            h=b(i,0);
            if(abs(g-h)<0.00001)
            {
                count++;
            }
            else
            {}
        }
        cout<<"-----------------"<<endl;
        cout<<count<<endl;
        cout<<"-------------------"<<endl;
        MatrixXd A1(count,n);
        MatrixXd A2((m-count), n);
        int count1=0,count2=0;
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
                count2++;
            }
        }
        cout<<"256 The Matrices are A1::\n"<<A1<<" A2::"<<endl<<A2<<endl;

        cout<<"============="<<count1<<"==========================="<<endl;

        if(count1>n ){
            episol=episol/2;
            b=degenerate(b_value,episol,m,n);
            cout<<" returned b="<<b<<endl;
            int dummy;
            //cin>>dummy;

            //getwchar();
            continue;
        }
        else if(count1<n){
            cout<<"Going into get_vetex\n"<<endl;
            x=get_vertex(A,b,x,m,n);
            b_value=b;
            //episol=episol/2;
            //b=degenerate(b_value,episol);
            continue;
        }
        cout<<"============================================="<<endl;
        MatrixXd B;
        A_copy=A1;
        B= -1*(A1.inverse());
        cout<<"================================"<<endl;
        cout<<B<<endl<<"Crosse degeneration diff\n"<<A1<<endl;
        cout<<"================================"<<endl;

        float k=0;
        int rows_a2[A2.rows()];
        float av[A2.rows()];
        for(int i=0;i<A2.rows();i++){
            rows_a2[i]=0;
            av[i]=0;
        }
        int z=0,stat=0;
        for(int i=0;i<B.cols();i++)
        {
            k =0;
            Bt=B.col(i);
            cout<<"==============================="<<endl;
            cout<<"Bt= \n"<<Bt<<endl;
            cout<<"+++++++++++++++++++++++++++++++"<<endl;
            if((ct*Bt)(0)>0)
            {
                cout<<" 302 You are not done\n"<<endl;
                stat++;
                B1=Bt;
                for(int j=0;j<A2.rows();j++)
                {
                    float sum=0;
                    for(int k=0;k<n;k++)
                    {
                        sum+=A2(j,k)*B1(k,0);
                    }

                    if(sum>0)
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
        if(stat==0){
            cout<<"Exiting from get_solution\n"<<endl;
            break;
        }
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
        }
        x=x+mini*B1;
        cout<<x<<endl;
        cout<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<mini<<"&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&&"<<endl;
    }

/*MatrixXd b_dash(A_copy.rows(),1);
int i,j;
for(i=0;i<A_copy.rows();i++)
    for(j=0;j<A.rows();j++)
        if(A_copy.row(i)==A.row(j))
            b_dash.row(i)=b.row(j);
x=A_copy.inverse()*b_dash;*/
cout<<"Final x=\n"<<x<<endl;
return x;
}

MatrixXd get_feasible_sol(MatrixXd A3, MatrixXd b3 ,int m, int n)
{
    int m1, n1;
    m1=m+1;n1=n+1;
    MatrixXd A(m1,n1),b(m1,1),c(n1,1),x(n1,1),Ax(m1,1),B1(n1,1),Bt(n1,1),ct(1,n1),sol(n1,1);

    for(int i=0;i<m1-1;i++){
        for(int j=0;j<n1;j++){
            if(j==n1-1)
                A(i,j)=1;
            if(j<n-1 && i==m-1)
                A(i,j)=0;
            if(j==n-1 && i==m1-1)
                A(i,j)=-1;
            A(i,j)=A3(i,j);
        }
    }
    for(int j=0;j<m1-1;j++)
        b(j,0)=b3(j,0);
    b(m1-1,0)= - b3.minCoeff();

    for(int i=0;i<n1-1;i++){
        x(i,0)=0;
        c(i,0)=0;
    }
    c(n1-1,0)=1;
    x(n1-1,0)=b3.minCoeff();
    bool rank_prb=false;
    MatrixXd cost_copy=c.transpose();
    FullPivLU<MatrixXd> lu_decomp(A);
    if(lu_decomp.rank()!=n1){
        MatrixXd dummy_a, dummy_b, dummy_c,dummy_x;
        dummy_a=A;dummy_b=b;dummy_c=c;dummy_x=x;
        n1=2*n1;m1=m1+n1;

        A.resize(m1,n1);b.resize(m1,1);c.resize(n1,1);x.resize(n1,1);
        Ax.resize(m1,1);ct.resize(1,n1);Bt.resize(n1,1);B1.resize(n1,1);
        A=convert_A(A,dummy_a,m1,n1);
        cout<<"Done A\n";
        b=convert_b(b,dummy_b,m1,n1);
        cout<<"Done b\n";
        c=convert_c(c,dummy_c,m1,n1);
        x=convert_x(x,dummy_x,m1,n1);

        rank_prb=true;
    }

    cout<<"A=\n"<<A<<endl;
    cout<<"b=\n"<<b<<endl;
    cout<<"c=\n"<<c<<endl;
    cout<<"x=\n"<<x<<endl;

    ct=c.transpose();
    cout<<"Going into get_Vertex"<<endl;

    x=get_vertex(A,b,x,m1,n1);
    int dummy;
    cout<<"x=\n"<<x<<endl;
    //cin>>dummy;


    x=get_solution(Ax,b,ct,A,x,Bt,B1,m1,n1);
    cout<<"Into the main\n";
    if(!rank_prb)
        sol=x;
    else{
        int i;
        {
            cout<<"INto the else\n";
            for(i=0;i<sol.rows();i++){

                sol(i,0)=(x(2*i,0)-x(2*i+1,0));
            }
        }

    }
    MatrixXd new_feasible(n,1);
    for(int i=0;i<n;i++)
        new_feasible(i,0)=sol(i,0);
    return new_feasible;
}





int main(int argc, char** argv) {
    int n,m;
    cout <<"Enter the no. of Rows and Columns of A:";
    cin>>m>>n;
    MatrixXd A(m,n),b(m,1),c(n,1),x(n,1),Ax(m,1),B1(n,1),Bt(n,1),ct(1,n),sol(n,1);
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

    char d;
    cout <<"Wanna enter initial feasible solution 'x':(y/N)"<<endl;
    cin>>d;
    if(d=='y'||d=='Y')
    {
        for(int l=0;l<n;l++)
            cin>>x(l,0);
    }
    else if(d=='N' || d=='n')
        x=get_feasible_sol(A,b,m,n);
    MatrixXd cost_copy=c.transpose();
    FullPivLU<MatrixXd> lu_decomp(A);
    if(lu_decomp.rank()!=n){
        MatrixXd dummy_a, dummy_b, dummy_c,dummy_x;
        dummy_a=A;dummy_b=b;dummy_c=c;dummy_x=x;
        n=2*n;m=m+n;

        A.resize(m,n);b.resize(m,1);c.resize(n,1);x.resize(n,1);
        Ax.resize(m,1);ct.resize(1,n);Bt.resize(n,1);B1.resize(n,1);
        A=convert_A(A,dummy_a,m,n);
        cout<<"Done A\n";
        b=convert_b(b,dummy_b,m,n);
        cout<<"Done b\n";
        c=convert_c(c,dummy_c,m,n);
        x=convert_x(x,dummy_x,m,n);

        rank_prb=true;
    }

    cout<<"A=\n"<<A<<endl;
    cout<<"b=\n"<<b<<endl;
    cout<<"c=\n"<<c<<endl;
    cout<<"x=\n"<<x<<endl;

    ct=c.transpose();
    cout<<"Going into get_Vertex"<<endl;

    x=get_vertex(A,b,x,m,n);
    int dummy;
    cout<<"x=\n"<<x<<endl;
    cin>>dummy;


    x=get_solution(Ax,b,ct,A,x,Bt,B1,m,n);
    cout<<"Into the main\n";
    if(!rank_prb)
        sol=x;
    else{
        int i;
        {
            cout<<"INto the else\n";
            for(i=0;i<sol.rows();i++){

                sol(i,0)=(x(2*i,0)-x(2*i+1,0));
            }
        }

    }
    cout<<"Optimum solution is "<<endl;
    cout <<sol<<endl;
    cout<<"Cost is:"<<cost_copy*sol<<endl;

    return 0;
}

