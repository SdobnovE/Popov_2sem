#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <stdlib.h>
using namespace std;

const double EPS = 1e-8;
string norm = "L2";
int m1, m2, n;
double mu, c;
const int T=1, X1=3, X2=3;
double tau, h1, h2;
double omega;

int ij2k (int i, int j);
void k2ij(int k, int& i, int& j);

double F1 (int i, int j, int t)
{
    return 0;
}
double F2 (int i, int j, int t)
{
    return 0;
}
double F3 (int i, int j, int t)
{
    return 0;
}

//ПРОВЕРИТЬ m1*m2 вместо K*K
/*
    Область
                 ___________
                |___|___|___|    
                |___|___|___|
     ___________|___|___|___|___________    m1,i
    |___|___|___|___|___|___|___|___|___|    
    |___|___|___|___|___|___|___|___|___|
    |___|___|___|___|___|___|___|___|___|
                    m2,j
*/



/*
    VEC = |____G____|____V1____|____V2_____| - массив (так как СЛУ одна)

*/

class Matrix
{
public:
    vector<vector<double>> Mat;// разреженная матрица(типа MSR, но другая)
    vector<vector<int>> Ind;//// индексы для разреженности
    vector<double> b;//вектор столбец для решения СЛУ
    int size;//размер матрицы (3 * K)


    Matrix();

    Matrix (int n): size(n)
    {
        Mat.resize(n);
        Ind.resize(n);
        b.resize(n);
    }
    
    

    void resize (int n)
    {
        size = n;
        Mat.resize(n);
        Ind.resize(n);
        b.resize(n);
    }

    vector<double> solveBICGS()
    {
        vector<double> d(size);
        vector<double> z_a(size);
        

        for (double& i: d)
            i = 1;
        

        for (double& i: z_a)
            i = 1; 

        
        

        vector<double> r(size), r_1(size);
        vector<double> temp(size),temp1(size);
        vector<double> p(size), u(size), q(size), s(size);


        


        sparse_matr_prod_vec(d, temp);
        for (int i = 0; i < size; i++)
        {
            u[i] = p[i] = r[i] = b[i] - temp[i];
        }

        int MAX_ITER = 100;
        int iter = 0;
        double alpha = 0, w = 0, beta = 0;
        double resid = -1;

        for (int I = 0; I < MAX_ITER; I++)
        {
            iter++;
            sparse_matr_prod_vec(p, temp);//temp = Ap_j
            alpha = scal_prod(r, z_a) / scal_prod(temp, z_a);

            for (int i = 0; i < size; i++)
                s[i] = r[i] - alpha * temp[i];
            
            
            sparse_matr_prod_vec(s, temp);//temp = As
            w = scal_prod(temp, s) / scal_prod(temp, temp);


            for (int i = 0; i < size; i++)
                d[i] = d[i] + alpha * p[i] + w * s[i];

            for (int i = 0; i < size; i++)
                r_1[i] = s[i] - w * temp[i];


            beta = scal_prod(r_1, z_a) / scal_prod(r, z_a) * alpha / w;


            sparse_matr_prod_vec(p, temp);//temp = Ap_j
            for (int i = 0; i < size; i++)
                p[i] = r_1[i] + beta * (p[i] - w * temp[i]);
            
            for (int i = 0; i < size; i++)
                r[i] = r_1[i];

            
            resid = sqrt(scal_prod(r,r));
            if (resid < EPS)
                break;
            
            //printf("\tresid:%e\n", resid);
            
        }

        sparse_matr_prod_vec (d, temp);

        
        for (int i = 0; i < size; i++)
            temp[i] -= b[i];
        double norm1 = sqrt (scal_prod(temp, temp));

        // printf ("\tResidual Solve Equations: %e\n", norm1);
        // printf ("\tITERATIONS: %d\n", iter);

        return d;
    }

    vector<double> solveEQU()//CGS
    {
        
        vector<double> d(size);
        vector<double> z_a(size);
        

        for (double& i: d)
            i = 1;
        

        for (double& i: z_a)
            i = 1; 

        
        

        vector<double> r(size), r_1(size);
        vector<double> temp(size),temp1(size);
        vector<double> p(size), u(size), q(size);



        


        sparse_matr_prod_vec(d, temp);
        for (int i = 0; i < size; i++)
        {
            u[i] = p[i] = r[i] = b[i] - temp[i];
        }

        int MAX_ITER = 2000;
        int iter = 0;
        double alpha = 0, beta = 0;
        double resid = -1;
        
        for (int tt = 0; tt < MAX_ITER ; tt++)
        {
            iter++;
            sparse_matr_prod_vec(p, temp);//temp = Ap_j
            alpha = scal_prod(r, z_a) / scal_prod(temp, z_a);
            
            for (int i = 0; i < size; i++)
                q[i] = u[i] - alpha * temp[i];

            for (int i = 0; i < size; i++)
                d[i] = d[i] + alpha * (u[i] + q[i]);
            
            
            for (int i = 0; i < size; i++)
                temp1[i] = u[i] + q[i];
            sparse_matr_prod_vec(temp1, temp);//A(u+q) = temp


            for (int i = 0; i < size; i++)
                r_1[i] = r[i] - alpha * temp[i];

            beta = scal_prod(r_1, z_a) / scal_prod(r, z_a);
            for (int i = 0; i < size; i++)
                u[i] = r_1[i] + beta * q[i];
            
            for (int i = 0; i < size; i++)
                p[i] = u[i] + beta * (q[i] + beta * p[i]);

            for (int i = 0; i < size; i++)
                r[i] = r_1[i];

            //printf("%e\n",sqrt(scal_prod(r,r)));
            if (resid  < 0)
            {
                
                resid = sqrt(scal_prod(r,r));
                //printf("\tresid:%e\n", resid);
                continue;
            }

            if (resid < EPS)
                break;
            resid = sqrt(scal_prod(r,r));
            //printf("\tresid:%e\n", resid);

        }

        
        sparse_matr_prod_vec (d, temp);

        
        for (int i = 0; i < size; i++)
            temp[i] -= b[i];
        double norm1 = sqrt (scal_prod(temp, temp));
        resid = sqrt(scal_prod(r,r));
        // printf ("\tresid: %e\n", resid);
        // printf ("\tResidual Solve Equations: %e\n", norm1);
        // printf ("\tITERATIONS: %d\n", iter);
        
        
        return d;
    }

    void fill_string(int type, int num_eq, int part, 
                    int i, int j, int t, vector<double> VEC1)
    {
        int K = ij2k((2*m1)/3, (2*m2)/3) + 1;
        Mat[num_eq].clear();
        Ind[num_eq].clear();
        // cout << num_eq % K << " " << part << " " << i << " " << j << " " << t << endl;
        double mu_wave = fabs(VEC1[0]);
        for (int i = 1; i < K; i++)
            mu_wave = fmax (fabs(VEC1[i]), mu_wave);
        mu_wave = exp(-mu_wave)*mu;

        if (part == 1)
        {
            switch (type)
            {
            case 0:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k(0, j) );
                //////////////////////////////////// G_0_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h1
                );
                Ind[num_eq].push_back(K + ij2k(1, j) );
                //////////////////////////////////// V1_1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h1
                );
                Ind[num_eq].push_back(K + ij2k(0, j) );
                //////////////////////////////////// V1_0_m2

                b[num_eq] = VEC1[ij2k(0, j)] * 1./tau + F1(0,j,t);
                break;
            
            case 1:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k(i, 0) );
                //////////////////////////////////// G_m1_0

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, 1) );
                //////////////////////////////////// V2_m1_1

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, 0) );
                //////////////////////////////////// V2_m1_0

                b[num_eq] = VEC1[ij2k(i, 0)] * 1./tau + F1(i,0,t);
                break;


            case 2:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, j) );
                //////////////////////////////////// V2_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1_M2-1

                b[num_eq] = VEC1[ij2k(i, j)] * 1./tau + F1(i,j,t);
                break;


            case 3:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k( i, j ) );
                //////////////////////////////////// G_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h1
                );
                Ind[num_eq].push_back(K + ij2k( i, j ) );
                //////////////////////////////////// V2_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h1
                );
                Ind[num_eq].push_back(K + ij2k( i-1, j ) );
                //////////////////////////////////// V2_m1_M2-1

                b[num_eq] = VEC1[ij2k(i, j)] * 1./tau + F1(i,j,t);
                break;

            case 4:
                 Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k( i, j ) );
                //////////////////////////////////// G_m1_0

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, 1 + j) );
                //////////////////////////////////// V2_m1_1

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k( i, j ) );
                //////////////////////////////////// V2_m1_0

                b[num_eq] = VEC1[ij2k(i, j)] * 1./tau + F1(i,j,t);
                


                break;

            case 5:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, j) );
                //////////////////////////////////// V2_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h2
                );
                Ind[num_eq].push_back(2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1_M2-1

                b[num_eq] = VEC1[ij2k(i, j)] * 1./tau+ F1(i,j,t);
                break;

            case 6:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                );
                Ind[num_eq].push_back( ij2k( i, j ) );
                //////////////////////////////////// G_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./h1
                );
                Ind[num_eq].push_back(K + ij2k( i, j ) );
                //////////////////////////////////// V2_m1_M2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./h1
                );
                Ind[num_eq].push_back(K + ij2k( i-1, j ) );
                //////////////////////////////////// V2_m1_M2-1

                b[num_eq] = VEC1[ij2k(i, j)] * 1./tau + F1(i,j,t);
                break;

            case 7:
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                // Mat[num_eq].push_back(1);
                // Ind[num_eq].push_back(num_eq);
                // b[num_eq] = log(rho(i,j,t));
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1.
                                        + tau * fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + tau * fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_m2


                                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -tau * ( fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                ////////////////////////////////////G_m1_m2-1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -tau*( fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2. * h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                ////////////////////////////////////G_m1-1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        tau*(VEC1[K + ij2k(i, j)] - fabs( VEC1[K + ij2k(i, j)] )) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                ////////////////////////////////////G_m1+1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        tau*(VEC1[2*K + ij2k(i, j)] - fabs( VEC1[2*K + ij2k(i, j)] )) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j+1) );
                ////////////////////////////////////G_m1_m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        tau/ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -tau/ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        tau/ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j+1) );
                //////////////////////////////////// V2_m1,m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -tau/ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1,m2-1

                b[num_eq] = VEC1[ij2k(i, j)] + tau*F1(i,j,t);
                break;

            }
        }
        else if (part == 2)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();
            if (type != 7)
            {//V1 = 0
                if (type == 6)
                {
                    Mat[num_eq].push_back(
                                            1.
                    );

                    Ind[num_eq].push_back( K + ij2k(i, j) );

                    b[num_eq] = omega;
                }
                else
                {
                    Mat[num_eq].push_back(
                                            1.
                    );

                    Ind[num_eq].push_back( K + ij2k(i, j) );

                    b[num_eq] = 0;
                }
            }
            else
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2./(h1 * h1) + 2./(h2 * h2))
                );
                Ind[num_eq].push_back( K + ij2k(i, j) );
                //////////////////////////////////// V2_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2) 
                );
                Ind[num_eq].push_back( K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave / (h1*h1)* 4./3
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V2_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave / (h1*h1)* 4./3
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V2_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j+1) );
                //////////////////////////////////// V2_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                //////////////////////////////////// G_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                //////////////////////////////////// G_m1_m2-1


                b[num_eq] = VEC1[K + ij2k(i, j)] / tau 
                            + F2(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        4./3. *( VEC1[2*K + ij2k(i-1, j)] - 2 * VEC1[2*K + ij2k(i, j)] + VEC1[2*K + ij2k(i+1, j)]) / (h1*h1)
                                        + ( VEC1[K + ij2k(i, j-1)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i, j+1)]) / (h2*h2)
                                      )
                            + mu/3. * exp( -VEC1[ij2k(i, j)] ) 
                                * ( VEC1[2*K + ij2k(i-1, j-1)] + VEC1[2*K + ij2k(i+1, j+1)] - VEC1[2*K + ij2k(i-1, j+1)] - VEC1[2*K + ij2k(i+1, j-1)] )
                                / (4 * h1 *h2)
                            ;
            }
            
        }
        else if (part == 3)
        {
            if (type != 7)
            {//V2 = 0
                if (type == 1)
                {
                    Mat[num_eq].push_back(
                                            1.
                    );

                    Ind[num_eq].push_back( 2*K + ij2k(i, j) );

                    b[num_eq] = omega;
                }
                else
                {
                    Mat[num_eq].push_back(
                                            1.
                    );

                    Ind[num_eq].push_back( 2*K + ij2k(i, j) );

                    b[num_eq] = 0;
                }
            }
            else
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2./(h2 * h2) + 2./(h1 * h1))
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j) );
                //////////////////////////////////// V2_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2) * 4./3.
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave / (h1*h1)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i-1, j) );
                //////////////////////////////////// V2_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave / (h1*h1)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i+1, j) );
                //////////////////////////////////// V2_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2) * 4./3.
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j+1) );
                //////////////////////////////////// V2_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c / (2*h2)
                );
                Ind[num_eq].push_back( ij2k(i, j+1) );
                //////////////////////////////////// G_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c / (2*h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                //////////////////////////////////// G_m1_m2-1


                b[num_eq] = VEC1[2*K + ij2k(i, j)] / tau 
                            + F3(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        ( VEC1[2*K + ij2k(i-1, j)] - 2 * VEC1[2*K + ij2k(i, j)] + VEC1[2*K + ij2k(i+1, j)]) / (h1*h1)
                                        + 4./3. *( VEC1[2*K + ij2k(i, j-1)] - 2 * VEC1[2*K + ij2k(i, j)] + VEC1[2*K + ij2k(i, j+1)]) / (h2*h2)
                                      )
                            + mu/3. * exp( -VEC1[ij2k(i, j)] ) 
                                * ( VEC1[K + ij2k(i-1, j-1)] + VEC1[K + ij2k(i+1, j+1)] - VEC1[K + ij2k(i-1, j+1)] - VEC1[K + ij2k(i+1, j-1)] )
                                / (4 * h1 *h2)
                            ;
            }   
        }
        
        

        
    }

    void fill_matrix_another()
    {
        int K = ij2k(m1/3*2, m2/3*2) + 1;
        vector<double> VEC1(3 * K);
        vector<double> VEC2(3 * K);

        for (int i = 0; i < K; i++)
        {
            int i1, j1;
            k2ij(i, i1, j1);
            
            VEC1[i] = 0;
            
            VEC1[i + K] = 0;
            VEC1[i + 2*K] = 0;
            
        }

        for(int i = 0; i <= m2/3; i++)//верхняя крышка
        {
            int i1,j1,k;
            i1 = 2*m1/3;
            j1 = i + m2/3;
            k = ij2k(i1, j1);
            VEC1[K + k] = omega;//u1

        }

        for(int i = 0; i <= m2/3; i++)//левая дальняя крышка
        {
            int i1,j1,k;
            i1 = i;
            j1 = 0;
            k = ij2k(i1, j1);
            VEC1[2*K + k] = omega;//u2

        }

        double resid_G=0, resid_V1=0, resid_V2=0;

        for (int t = 1; t <= n; t++)
        {
            

            for (int k = 0; k < K; k++)
            {
                int i1, j1;
                k2ij(k, i1, j1);
                int type = 0;
                
                if (i1 == 0)
                    type = 0;//нижняя крышка;
                else if (j1 == 0)
                    type = 1; //левая дальняя крышка
                else if (j1 == m2)
                    type = 2; //правая дальняя крышка
                else if (i1 == m1/3 && (j1 <= m2/3 || j1 >= 2*m2/3))
                    type = 3;//нижние верхние крышки
                else if (i1 > m1/3 && j1 == m2/3)
                    type = 4; // левая ближняя крышка
                else if (i1 > m1/3 && j1 == 2*m2/3)
                    type = 5; // правая ближжняя крышка
                else if (i1 == 2*m1/3)
                    type = 6;//верхняя крышка
                else
                    type = 7;//внутренняя точка
                
                

                fill_string(type, k, 1, i1, j1, t, VEC1);
                
                fill_string(type, K + k, 2, i1, j1, t, VEC1);

                
                fill_string(type, 2*K + k, 3, i1, j1, t, VEC1);

                
            }

            

            auto res = solveEQU();
            
            
            for (int i = 0; i < res.size(); i++)
                VEC1[i] = res[i];

            

        }


        


        vector<vector<double>> RESULT(m1+1), INIT(m1+1);
        
        for (int i = 0; i < m1+1; i++)
        {
            RESULT[i].resize(m2+1);
            INIT[i].resize(m2+1);
            for (int j = 0; j < m2+1; j++)
            {
                RESULT[i][j] = 0.;
                INIT[i][j] = 0.;
            }
        }


        for (int i = 0; i < K; i++)
        {
            
            int i1, j1;
            k2ij(i, i1, j1);
            
            RESULT[i1][j1] = VEC1[i];
            
        }
        
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < m2; j++)
                cout << RESULT[i][j] << " ";
            cout << endl;
        }
        















        for (int i = 0; i < m1+1; i++)
        {
            RESULT[i].resize(m2+1);
            INIT[i].resize(m2+1);
            for (int j = 0; j < m2+1; j++)
            {
                RESULT[i][j] = 0.;
                INIT[i][j] = 0.;
            }
        }


        for (int i = 0; i < K; i++)
        {
            
            int i1, j1;
            k2ij(i, i1, j1);
            
            RESULT[i1][j1] = VEC1[i+K];
            
        }
        
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < m2; j++)
                cout << RESULT[i][j] << " ";
            cout << endl;
        }
        












        for (int i = 0; i < m1+1; i++)
        {
            RESULT[i].resize(m2+1);
            INIT[i].resize(m2+1);
            for (int j = 0; j < m2+1; j++)
            {
                RESULT[i][j] = 0.;
                INIT[i][j] = 0.;
            }
        }


        for (int i = 0; i < K; i++)
        {
            
            int i1, j1;
            k2ij(i, i1, j1);
            
            RESULT[i1][j1] = VEC1[i+2*K];
            
        }
        
        for (int i = 0; i < m1; i++)
        {
            for (int j = 0; j < m2; j++)
                cout << RESULT[i][j] << " ";
            cout << endl;
        }
        
        
    }


    
    void count_residual(vector<double> VEC1, vector<double> VEC2, double& resid_G, double& resid_V1, double& resid_V2)
    {
        int K = ij2k(m1/3*2, m2/3*2) + 1;

        if (norm == "C")
            {
                resid_V2 = resid_V1 = resid_G = 0;
                for (int i = 0; i < K; i++)
                {
                    if (fabs (VEC1[i] - VEC2[i]) > resid_G)
                        resid_G = fabs (VEC2[i] - VEC1[i]);
                }
                
                for (int i = K; i < 2*K; i++)
                {
                    if (fabs (VEC1[i] - VEC2[i]) > resid_V1)
                        resid_V1 = fabs (VEC2[i] - VEC1[i]);
                }

                for (int i = 2*K; i < 3*K; i++)
                {
                    if (fabs (VEC1[i] - VEC2[i]) > resid_V2)
                        resid_V2 = fabs (VEC2[i] - VEC1[i]);
                }
            }

            if (norm == "L2")
            {
                resid_V2 = resid_V1 = resid_G = 0;
                for(int k = 0; k < K; k++)
                {
                    int i1, j1;
                    k2ij(k, i1, j1);
                    int type = 0;
                    
                    if (i1 == 0)
                        type = 0;//нижняя крышка;
                    else if (j1 == 0)
                        type = 1; //левая дальняя крышка
                    else if (j1 == m2)
                        type = 2; //правая дальняя крышка
                    else if (i1 == m1/3 && (j1 <= m2/3 || j1 >= 2*m2/3))
                        type = 3;//нижние верхние крышки
                    else if (i1 > m1/3 && j1 == m2/3)
                        type = 4; // левая ближняя крышка
                    else if (i1 > m1/3 && j1 == 2*m2/3)
                        type = 5; // правая ближжняя крышка
                    else if (i1 == 2*m1/3)
                        type = 6;//верхняя крышка
                    else
                        type = 7;//внутренняя точка
                    
                    if (type == 7)
                    {
                        resid_G += (VEC1[k] - VEC2[k]) * (VEC1[k] - VEC2[k]);
                        resid_V1 += (VEC1[K+k] - VEC2[K+k]) * (VEC1[K+k] - VEC2[K+k]);
                        resid_V2 += (VEC1[2*K+k] - VEC2[2*K+k]) * (VEC1[2*K+k] - VEC2[2*K+k]);
                    }
                    else
                    {
                        resid_G += 0.5*(VEC1[k] - VEC2[k]) * (VEC1[k] - VEC2[k]);
                        resid_V1 += 0.5*(VEC1[K+k] - VEC2[K+k]) * (VEC1[K+k] - VEC2[K+k]);
                        resid_V2 += 0.5*(VEC1[2*K+k] - VEC2[2*K+k]) * (VEC1[2*K+k] - VEC2[2*K+k]);
                    }

                    
                }
                
                resid_G = sqrt(resid_G*h1*h2);
                resid_V1 = sqrt(resid_V1*h1*h2);
                resid_V2 = sqrt(resid_V2*h1*h2);

                
            }


            if (norm == "W")
            {
                resid_V2 = resid_V1 = resid_G = 0;
                double t_G = 0., t_V1 = 0, t_V2 = 0;

                for(int k = 0; k < K; k++)
                {
                    int i1, j1;
                    k2ij(k, i1, j1);
                    int type = 0;
                    
                    if (i1 == 0)
                        type = 0;//нижняя крышка;
                    else if (j1 == 0)
                        type = 1; //левая дальняя крышка
                    else if (j1 == m2)
                        type = 2; //правая дальняя крышка
                    else if (i1 == m1/3 && (j1 <= m2/3 || j1 >= 2*m2/3))
                        type = 3;//нижние верхние крышки
                    else if (i1 > m1/3 && j1 == m2/3)
                        type = 4; // левая ближняя крышка
                    else if (i1 > m1/3 && j1 == 2*m2/3)
                        type = 5; // правая ближжняя крышка
                    else if (i1 == 2*m1/3)
                        type = 6;//верхняя крышка
                    else
                        type = 7;//внутренняя точка
                    if (type == 7)
                    {
                        resid_G += (VEC1[k] - VEC2[k]) * (VEC1[k] - VEC2[k]);
                        resid_V1 += (VEC1[K+k] - VEC2[K+k]) * (VEC1[K+k] - VEC2[K+k]);
                        resid_V2 += (VEC1[2*K+k] - VEC2[2*K+k]) * (VEC1[2*K+k] - VEC2[2*K+k]);
                    }
                    else
                    {
                        resid_G += 0.5*(VEC1[k] - VEC2[k]) * (VEC1[k] - VEC2[k]);
                        resid_V1 += 0.5*(VEC1[K+k] - VEC2[K+k]) * (VEC1[K+k] - VEC2[K+k]);
                        resid_V2 += 0.5*(VEC1[2*K+k] - VEC2[2*K+k]) * (VEC1[2*K+k] - VEC2[2*K+k]);
                    }


                    
                    
                    for (int i = 0; i <= m1/3; i++)// (res - true)_x2
                    {
                        for (int j = 0; j < m2; j++)
                        {
                            int k1 = ij2k(i, j);
                            int k2 = ij2k(i, j+1);
                            
                            t_G += ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h2
                                * ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h2;
                            t_V1 += ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h2
                                * ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h2;
                            t_V2 += ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h2
                                * ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h2;
                        }
                    }

                    for (int i = m1/3 + 1; i <= 2*m1/3; i++)// (res - true)_x2
                    {
                        for (int j = m2/3; j < 2*m2/3; j++)
                        {
                            int k1 = ij2k(i, j);
                            int k2 = ij2k(i, j+1);
                            
                            t_G += ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h2
                                * ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h2;
                            t_V1 += ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h2
                                * ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h2;
                            t_V2 += ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h2
                                * ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h2;
                        }
                    }

                    for (int i = 0; i < m1/3; i++)// (res - true)_x2
                    {
                        for (int j = 0; j <= m2; j++)
                        {
                            int k1 = ij2k(i, j);
                            int k2 = ij2k(i+1, j);
                            
                            t_G += ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h1
                                * ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h1;
                            t_V1 += ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h1 
                                * ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h1;
                            t_V2 += ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h1
                                * ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h1;
                        }
                    }

                    for (int i = m1/3; i < 2*m1/3; i++)// (res - true)_x2
                    {
                        for (int j = m2/3; j <= 2*m2/3; j++)
                        {
                            int k1 = ij2k(i, j);
                            int k2 = ij2k(i+1, j);
                            
                            t_G += ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h1
                                * ((VEC1[k1] - VEC2[k1]) - (VEC1[k2] - VEC2[k2])) / h1;
                            t_V1 += ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h1 
                                * ((VEC1[K+k1] - VEC2[K+k1]) - (VEC1[K+k2] - VEC2[K+k2])) / h1;
                            t_V2 += ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h1
                                * ((VEC1[2*K+k1] - VEC2[2*K+k1]) - (VEC1[2*K+k2] - VEC2[2*K+k2])) / h1;
                        }
                    }


                    
                }
                

                resid_G += t_G;
                resid_V1 += t_V1;
                resid_V2 += t_V2;

                resid_G = sqrt(resid_G*h1*h2);
                resid_V1 = sqrt(resid_V1*h1*h2);
                resid_V2 = sqrt(resid_V2*h1*h2);
                

                
            }
    }

    

private:
    void sparse_matr_prod_vec (const vector<double>& x, vector<double>& res)//Ax=res
    {
        res.clear();
        res.resize(Mat.size());
        
        for (double& i : res)
            i = 0.;
        
        for (int i = 0; i < size; i++)
        {
            for (int j = 0; j < Mat[i].size(); j++)
            {
                res[i] += Mat[i][j] * x[Ind[i][j]];
            }
        }
    }

    double scal_prod (const vector<double>& vec1, const vector<double>& vec2)
    {
        if (vec1.size() != vec2.size())
            throw runtime_error("different dimentions in scal_prod");

        double res = 0;
        for (int i = 0; i < vec1.size(); i++)
        {
            res += vec1[i] * vec2[i];
        }

        return res;
    }

};
Matrix Mat(1);






int ij2k (int i, int j)
{
    if (i <= m1/3)
    {
        return (m2 + 1) * i + j;
    }
    else
    {   
        return (m2 + 1) * (m1 / 3 + 1) + (i - m1/3-1) * (m2/3 + 1) + (j - m2/3);
    }
    
}


void k2ij(int k, int& i, int& j)
{
    if (k < (m2 + 1) * (m1 / 3 + 1))
    {
        i = k / (m2 + 1);
        j = k % (m2 + 1);
    }
    else
    {
        i = m1/3 + 1;
        j = m2/3;
        k -= (m2 + 1) * (m1 / 3 + 1);

        i += k / (m2/3 + 1);
        j += k % (m2/3 + 1);
    }
    
}


int main(int argc, char *argv[])
{
    if (argc != 8)
    {
        cout << "USAGE: ./a.out m1 m2 n mu c norm omega\n";
        return -1;
    }

    m1 = atoi(argv[1]);
    m2 = atoi(argv[2]);
    n = atoi(argv[3]);
    mu = atof(argv[4]);
    c = atof(argv[5]);
    norm = argv[6];
    omega = atof(argv[7]);
    Mat.resize(3*ij2k(2*m1/3, 2*m2/3) + 3);

    tau = 0.01;
    h1 = static_cast<double>(X1) / m1;
    h2 = static_cast<double>(X2) / m2;

    
    Mat.fill_matrix_another();

    return 0;
}
