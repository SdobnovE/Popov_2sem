#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <stdlib.h>
using namespace std;

const double EPS = 1e-16;

int m1, m2, n;
double mu, c;
const int T=1, X1=1, X2=1;
double tau, h1, h2;

int ij2k (int i, int j);
void k2ij(int k, int& i, int& j);
double u1(int i, int j, int t);
double u2(int i, int j, int t);
double rho(int i, int j, int t);
double F1 (int i, int j, int t);
double F2 (int i, int j, int t);
double F3 (int i, int j, int t);
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
    
    void pirnt_normal()
    {
        FILE* f;
        f = fopen("file", "w");

        vector<double> res(size);
        for (int i = 0; i < size; i++)
        {
            res.clear();
            res.resize(size);
            for (int j = 0; j < Mat[i].size(); j++)
            {
                res[Ind[i][j]] = Mat[i][j];
            }
            for (auto i : res)
                fprintf(f, "%e ", i);
            fprintf(f, "\n");
        }
        fclose(f);
        f = fopen("b", "w");
        for (auto i : b)
                fprintf(f, "%e\n", i);
        fclose(f);
    }

    void resize (int n)
    {
        size = n;
        Mat.resize(n);
        Ind.resize(n);
        b.resize(n);
    }

    vector<double> solveBICGS()//CGS
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
            
            printf("\tresid:%e\n", resid);
            
        }

        sparse_matr_prod_vec (d, temp);

        
        for (int i = 0; i < size; i++)
            temp[i] -= b[i];
        double norm1 = sqrt (scal_prod(temp, temp));

        printf ("\tResidual Solve Equations: %e\n", norm1);
        printf ("\tITERATIONS: %d\n", iter);

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

        int MAX_ITER = 100;
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
                printf("\tresid:%e\n", resid);
                continue;
            }

            if (resid < EPS)
                break;
            resid = sqrt(scal_prod(r,r));
            printf("\tresid:%e\n", resid);

        }

        
        sparse_matr_prod_vec (d, temp);

        
        for (int i = 0; i < size; i++)
            temp[i] -= b[i];
        double norm1 = sqrt (scal_prod(temp, temp));

        printf ("\tResidual Solve Equations: %e\n", norm1);
        printf ("\tITERATIONS: %d\n", iter);
        
        
        return d;
    }

    void fill_string(int type, int num_eq, int part)
    {
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
            
            default:
                break;
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
            
            VEC1[i] = log(rho(i1, j1, 0));
            // cout << i1 << " " << j1 << log(rho(i1, j1, 0)) << endl;
            VEC1[i + K] = u1(i1, j1, 0);
            VEC1[i + 2*K] = u2(i1, j1, 0);
            
        }

        // for (auto i: VEC1)
        //     cout << "\t" << i << endl;

        for (int i = 0; i < K; i++)
        {
            int i1, j1;
            k2ij(i, i1, j1);
            VEC2[i] = log(rho(i1, j1, 1));
            VEC2[i + K] = u1(i1, j1, 1);
            VEC2[i + 2*K] = u2(i1, j1, 1);
            
        }

        int t = 1;

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

            fill_string(type, k, 1);
        }

    }
    


    void fill_matrix()
    {
        int K = ij2k(m1/3*2, m2/3*2) + 1;
        vector<double> VEC1(3 * K);
        vector<double> VEC2(3 * K);

        for (int i = 0; i < K; i++)
        {
            int i1, j1;
            k2ij(i, i1, j1);
            
            VEC1[i] = log(rho(i1, j1, 0));
            // cout << i1 << " " << j1 << log(rho(i1, j1, 0)) << endl;
            VEC1[i + K] = u1(i1, j1, 0);
            VEC1[i + 2*K] = u2(i1, j1, 0);
            
        }

        // for (auto i: VEC1)
        //     cout << "\t" << i << endl;

        for (int i = 0; i < K; i++)
        {
            int i1, j1;
            k2ij(i, i1, j1);
            VEC2[i] = log(rho(i1, j1, 1));
            VEC2[i + K] = u1(i1, j1, 1);
            VEC2[i + 2*K] = u2(i1, j1, 1);
            
        }

        int t = 1;
        int num_eq = 0; 
        
        double mu_wave = fabs(VEC1[0]);
        for (int i = 1; i < K; i++)
            mu_wave = max (fabs(VEC1[i]), mu_wave);
        mu_wave = exp(-mu_wave)*mu; 
        
        ///////////////////////////////////////////////////////////////////////////////////////////////// 1-ое уравнение
        for (int i = 1; i < m1/3; i++)
        {
            for (int j = 1; j < m2; j++)
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_m2


                                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                ////////////////////////////////////G_m1_m2-1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                ////////////////////////////////////G_m1-1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[K + ij2k(i, j)] - fabs( VEC1[K + ij2k(i, j)] )) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                ////////////////////////////////////G_m1+1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[2*K + ij2k(i, j)] - fabs( VEC1[2*K + ij2k(i, j)] )) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j+1) );
                ////////////////////////////////////G_m1_m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j+1) );
                //////////////////////////////////// V2_m1,m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1,m2-1

                b[num_eq] = VEC1[ij2k(i, j)] / tau + F1(i,j,t);
                num_eq++;


            }
        }

        for (int i = 0; i < m1/3; i++)
        {
            for (int j = 1; j < m2/3; j++)
            {
                i += m1/3;
                j += m2/3;

                Mat[num_eq].clear();
                Ind[num_eq].clear();

                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_m2


                                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                ////////////////////////////////////G_m1_m2-1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                ////////////////////////////////////G_m1-1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[K + ij2k(i, j)] - fabs( VEC1[K + ij2k(i, j)] )) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                ////////////////////////////////////G_m1+1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[2*K + ij2k(i, j)] - fabs( VEC1[2*K + ij2k(i, j)] )) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j+1) );
                ////////////////////////////////////G_m1_m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./ (2*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1,m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j+1) );
                //////////////////////////////////// V2_m1,m2+1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -1./ (2*h2)
                );
                Ind[num_eq].push_back( 2*K + ij2k(i, j-1) );
                //////////////////////////////////// V2_m1,m2-1

                b[num_eq] = VEC1[ij2k(i, j)] / tau + F1(i,j,t);

                num_eq++;
                i -= m1/3;
                j -= m2/3;
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


        
        


        ///////////////////////////////////////////////////////////////////////////////////////////////// 3-е уравнение
        for (int i = 1; i < m1/3; i++)
        {
            for (int j = 1; j < m2; j++)
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h1 * h1) + 2/(h2 * h2))
                );
                Ind[num_eq].push_back( K + ij2k(i, j) );
                //////////////////////////////////// V1_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j-1) );
                //////////////////////////////////// V1_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j+1) );
                //////////////////////////////////// V1_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                //////////////////////////////////// G_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                //////////////////////////////////// G_m1-1_m2


                b[num_eq] = VEC1[K + ij2k(i, j)] / tau 
                            + F2(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        4./3. * ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
                                        + ( VEC1[K + ij2k(i, j-1)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i, j+1)]) / (h2*h2)
                                      )
                            + mu/3. * exp( -VEC1[ij2k(i, j)] ) 
                                * ( VEC1[2*K + ij2k(i-1, j-1)] + VEC1[2*K + ij2k(i+1, j+1)] - VEC1[2*K + ij2k(i-1, j+1)] - VEC1[2*K + ij2k(i+1, j-1)] )
                                / (4 * h1 *h2)
                            ;

                num_eq++;


            }
        }

        for (int i = 0; i < m1/3; i++)
        {
            for (int j = 1; j < m2/3; j++)
            {
                i += m1/3;
                j += m2/3;

                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h1 * h1) + 2/(h2 * h2))
                );
                Ind[num_eq].push_back( K + ij2k(i, j) );
                //////////////////////////////////// V1_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j-1) );
                //////////////////////////////////// V1_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-fabs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j+1) );
                //////////////////////////////////// V1_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                //////////////////////////////////// G_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c / (2*h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                //////////////////////////////////// G_m1-1_m2


                b[num_eq] = VEC1[K + ij2k(i, j)] / tau 
                            + F2(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        4./3. * ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
                                        + ( VEC1[K + ij2k(i, j-1)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i, j+1)]) / (h2*h2)
                                      )
                            + mu/3. * exp( -VEC1[ij2k(i, j)] ) 
                                * ( VEC1[2*K + ij2k(i-1, j-1)] + VEC1[2*K + ij2k(i+1, j+1)] - VEC1[2*K + ij2k(i-1, j+1)] - VEC1[2*K + ij2k(i+1, j-1)] )
                                / (4 * h1 *h2)
                            ;

                num_eq++;
                i -= m1/3;
                j -= m2/3;
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////
        
        
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////4-ое уравнение
        for (int i = 1; i < m1/3; i++)
        {
            for (int j = 1; j < m2; j++)
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h2 * h2) + 2/(h1 * h1))
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

                num_eq++;


            }
        }  


        for (int i = 0; i < m1/3; i++)
        {
            for (int j = 1; j < m2/3; j++)
            {
                i += m1/3;
                j += m2/3;

                Mat[num_eq].clear();
                Ind[num_eq].clear();

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        1./tau
                                        + fabs( VEC1[K + ij2k(i, j)] ) / h1
                                        + fabs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h2 * h2) + 2/(h1 * h1))
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

                num_eq++;
                i -= m1/3;
                j -= m2/3;
            }
        }  
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////

        



        
        
        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// уравнение V1=0 V2=0 на границе

        for (int j = 0; j <= m2; j++)//нижняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back(K + j);

            b[num_eq] = 0;
            num_eq++;
            
        }

        for (int i = 1; i <= m1/3; i++)//левая дальняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(i, 0) );

            b[num_eq] = 0;
            num_eq++;
        }


        for (int i = 1; i <= m1/3; i++)//правая дальняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(i, m2) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int i = 0; i <= m1/3; i++)//левая ближняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(i + m1/3, m2/3) );

            b[num_eq] = 0;
            num_eq++;
        }


        for (int i = 0; i <= m1/3; i++)//правая ближняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(i + m1/3, 2*m2/3) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя центральная часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(2 *m1/3, m2/3 + j) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя левая часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(m1/3, j) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя правая часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( K + ij2k(m1/3, 2*m2/3 + j) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 0; j <= m2; j++)//нижняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back(2*K + j);

            b[num_eq] = 0;
            num_eq++;
            
        }

        for (int i = 1; i <= m1/3; i++)//левая дальняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(i, 0) );

            b[num_eq] = 0;
            num_eq++;
        }


        for (int i = 1; i <= m1/3; i++)//правая дальняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(i, m2) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int i = 0; i <= m1/3; i++)//левая ближняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(i + m1/3, m2/3) );

            b[num_eq] = 0;
            num_eq++;
        }


        for (int i = 0; i <= m1/3; i++)//правая ближняя часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(i + m1/3, 2*m2/3) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя центральная часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(2 *m1/3, m2/3 + j) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя левая часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(m1/3, j) );

            b[num_eq] = 0;
            num_eq++;
        }

        for (int j = 1; j < m2/3; j++)//верхняя правая часть
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            Mat[num_eq].push_back(1);
            Ind[num_eq].push_back( 2*K + ij2k(m1/3, 2*m2/3 + j) );

            b[num_eq] = 0;
            num_eq++;
        }

        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////


        /////////////////////////////////////////////////////////////////////////////////////////////////////////////// 2-ое уравнение

        for (int j = 0; j <= m2; j++)//нижняя часть(1-ое уравнение)
        {
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
            

            num_eq++;

        }

        for (int i = 1; i <= m1/3; i++)//левая дальняя часть(2-ое уравнение)
        {
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
            

            num_eq++;

        }


        for (int i = 0; i <= m1/3; i++)//левая ближняя часть(2-ое уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(i + m1/3, m2/3) );
            //////////////////////////////////// G_m1_0

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i + m1/3, 1 + m2/3) );
            //////////////////////////////////// V2_m1_1

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i + m1/3, m2/3) );
            //////////////////////////////////// V2_m1_0

            b[num_eq] = VEC1[ij2k(i + m1/3, m2/3)] * 1./tau + F1(i + m1/3,m2/3,t);
            

            num_eq++;

        }

        for (int i = 1; i <= m1/3; i++)//правая дальняя часть(4-ое уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(i, m2) );
            //////////////////////////////////// G_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i, m2) );
            //////////////////////////////////// V2_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i, m2-1) );
            //////////////////////////////////// V2_m1_M2-1

            b[num_eq] = VEC1[ij2k(i, m2)] * 1./tau+ F1(i,m2,t);
            

            num_eq++;

        }


        for (int i = 0; i <= m1/3; i++)//правая ближняя часть(4-ое уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(i + m1/3, 2*m2/3) );
            //////////////////////////////////// G_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i + m1/3, 2*m2/3) );
            //////////////////////////////////// V2_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h2
            );
            Ind[num_eq].push_back(2*K + ij2k(i + m1/3, 2*m2/3-1) );
            //////////////////////////////////// V2_m1_M2-1

            b[num_eq] = VEC1[ij2k(i + m1/3, 2*m2/3)] * 1./tau+ F1(i + m1/3,2*m2/3,t);
            

            num_eq++;

        }

        for (int j = 1; j < m2/3; j++)//верхняя центральная часть(3-е уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(2*m1/3, j + m2/3) );
            //////////////////////////////////// G_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h1
            );
            Ind[num_eq].push_back(K + ij2k(2*m1/3, j + m2/3) );
            //////////////////////////////////// V2_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h1
            );
            Ind[num_eq].push_back(K + ij2k(2*m1/3-1, j + m2/3) );
            //////////////////////////////////// V2_m1_M2-1

            b[num_eq] = VEC1[ij2k(2*m1/3, j + m2/3)] * 1./tau + F1(2*m1/3,j + m2/3,t);
            

            num_eq++;
   
        }

        for (int j = 1; j < m2/3; j++)//верхняя левая часть(3-е уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(m1/3, j) );
            //////////////////////////////////// G_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h1
            );
            Ind[num_eq].push_back(K + ij2k(m1/3, j) );
            //////////////////////////////////// V2_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h1
            );
            Ind[num_eq].push_back(K + ij2k(m1/3-1, j) );
            //////////////////////////////////// V2_m1_M2-1

            b[num_eq] = VEC1[ij2k(m1/3, j)] * 1./tau+ F1(m1/3,j,t);
            

            num_eq++;
   
        }

        for (int j = 1; j < m2/3; j++)//верхняя правая часть(3-е уравнение)
        {
            Mat[num_eq].clear();
            Ind[num_eq].clear();

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./tau
            );
            Ind[num_eq].push_back( ij2k(m1/3, j + 2*m2/3) );
            //////////////////////////////////// G_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    1./h1
            );
            Ind[num_eq].push_back(K + ij2k(m1/3, j + 2*m2/3) );
            //////////////////////////////////// V2_m1_M2

            ////////////////////////////////////
            Mat[num_eq].push_back(
                                    -1./h1
            );
            Ind[num_eq].push_back(K + ij2k(m1/3-1, j + 2*m2/3) );
            //////////////////////////////////// V2_m1_M2-1

            b[num_eq] = VEC1[ij2k(m1/3, j + 2*m2/3)] * 1./tau + F1(m1/3,j + 2*m2/3,t);
            

            num_eq++;
   
        }



        //this->pirnt_normal();
        // for (auto& i : Mat)
        // {
        //     for (auto& j : i)
        //         cout << j <<"\t";
        //     cout << endl;
        // }

        //this->pirnt_normal();
        this->solveBICGS();

        
        
    //     sparse_matr_prod_vec(VEC2, VEC1);
    //     VEC1 = {1.87509250e+00,  1.91467029e+00,  1.94694961e+00,  1.97774237e+00,
    //     2.00213649e+00,  2.02862200e+00,  2.05429365e+00,  2.07875702e+00,
    //     2.10418960e+00,  2.12517673e+00,  2.14635476e+00,  2.17090162e+00,
    //     2.22626750e+00,  1.87833705e+00,  1.91725937e+00,  1.95169007e+00,
    //     1.98268942e+00,  2.01225793e+00,  2.03535521e+00,  2.06177039e+00,
    //     2.08577710e+00,  2.11482552e+00,  2.13930270e+00,  2.16765515e+00,
    //     2.20496066e+00,  2.27638312e+00,  1.88274505e+00,  1.92317983e+00,
    //     1.95697547e+00,  1.99033683e+00,  2.00690737e+00,  2.03515570e+00,
    //     2.06415827e+00,  2.09326713e+00,  2.12938867e+00,  2.15515388e+00,
    //     2.18621943e+00,  2.23160332e+00,  2.33073490e+00,  1.88836223e+00,
    //     1.92939051e+00,  1.97151460e+00,  2.00754228e+00,  2.04555243e+00,
    //     2.04305146e+00,  2.07518878e+00,  2.10485530e+00,  2.17921745e+00,
    //     2.21381458e+00,  2.25181230e+00,  2.32201411e+00,  2.39301903e+00,
    //     1.90828417e+00,  1.95157037e+00,  2.00418645e+00,  2.05569135e+00,
    //     1.96687608e+00,  2.01482481e+00,  2.07520719e+00,  2.14612575e+00,
    //     2.28436346e+00,  2.30732033e+00,  2.34164792e+00,  2.37498231e+00,
    //     2.26421000e+00,  1.95345657e+00,  2.00882608e+00,  2.07826805e+00,
    //     2.16423074e+00,  2.30848244e+00,  1.93744872e+00,  2.01028980e+00,
    //     2.08747468e+00,  2.18541318e+00,  2.33279276e+00,  1.93087987e+00,
    //     2.04608402e+00,  2.14002802e+00,  2.27039221e+00,  2.35919442e+00,
    //     2.03450167e+00,  2.19520287e+00,  2.28231221e+00,  2.35737010e+00,
    //     2.15284000e+00,  0.00000000e+00, -2.39610989e-15,  2.95326839e-16,
    //    -2.69654693e-15,  1.83456569e-16, -1.72248984e-15,  2.42165226e-15,
    //    -8.58917598e-16,  2.81576804e-16, -1.70219799e-15,  1.17551526e-16,
    //    -1.66383984e-15,  0.00000000e+00,  0.00000000e+00,  1.19804196e-03,
    //     7.83038931e-03,  1.39792977e-02,  2.44726733e-02,  3.07546703e-02,
    //     3.56788538e-02,  3.95988088e-02,  4.03078954e-02,  4.31966011e-02,
    //     4.36135764e-02,  3.83675440e-02,  0.00000000e+00,  0.00000000e+00,
    //     4.47776803e-03,  1.56368663e-02,  2.77965916e-02,  3.93326787e-02,
    //     6.34787515e-02,  7.44554002e-02,  8.36258807e-02,  7.58578339e-02,
    //     8.06354347e-02,  8.56624626e-02,  7.98000925e-02,  0.00000000e+00,
    //     0.00000000e+00,  1.86869862e-03,  1.49864489e-02,  2.90121793e-02,
    //     9.23534497e-02,  1.10245984e-01,  1.16834994e-01,  1.17189072e-01,
    //     9.04990829e-02,  1.01677000e-01,  1.14154591e-01,  1.27964814e-01,
    //     0.00000000e+00,  0.00000000e+00,  9.47885521e-17, -3.36795941e-16,
    //    -1.27990385e-16,  2.94680255e-17,  1.36632302e-01,  1.59068344e-01,
    //     1.58982648e-01,  5.36792390e-16, -3.13722442e-16, -3.79772702e-16,
    //     1.27941086e-15,  0.00000000e+00,  5.29800467e-16,  1.51456406e-01,
    //     1.92300494e-01,  2.08980261e-01,  0.00000000e+00, -5.40309920e-16,
    //     1.58969624e-01,  2.07676742e-01,  2.27676185e-01,  0.00000000e+00,
    //    -1.24954709e-16,  1.27803703e-01,  1.84198873e-01,  2.30766770e-01,
    //     0.00000000e+00,  0.00000000e+00,  1.37256873e-15,  1.10923881e-16,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00,  1.78756969e-16,
    //     9.47128656e-03,  1.80703382e-02,  2.52702946e-02,  3.35895514e-02,
    //     3.95755990e-02,  4.35697777e-02,  4.65416511e-02,  5.02476008e-02,
    //     5.41704075e-02,  5.18163459e-02,  3.56464492e-02,  0.00000000e+00,
    //    -7.72479569e-17,  1.50791194e-02,  2.97472888e-02,  4.84865361e-02,
    //     6.19429944e-02,  6.72723114e-02,  7.62052045e-02,  8.23650478e-02,
    //     9.22037573e-02,  1.02527280e-01,  1.01670195e-01,  7.85574010e-02,
    //     0.00000000e+00, -4.41853297e-16,  1.64052664e-02,  3.59680108e-02,
    //     5.54434913e-02,  8.80764562e-02,  1.07816751e-01,  1.10617005e-01,
    //     1.07364840e-01,  1.16491158e-01,  1.41830918e-01,  1.41351539e-01,
    //     1.32877367e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00, -2.47174413e-16,  9.51805824e-02,
    //     1.14454071e-01,  1.02885127e-01,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00, -4.43058203e-16,
    //     1.08383431e-01,  1.33946936e-01,  1.27212438e-01,  0.00000000e+00,
    //    -4.93331004e-17,  1.19877948e-01,  1.55222310e-01,  1.56208590e-01,
    //     0.00000000e+00, -5.36548716e-17,  1.17394297e-01,  1.62147671e-01,
    //     1.92026916e-01,  0.00000000e+00,  0.00000000e+00,  0.00000000e+00,
    //     0.00000000e+00,  0.00000000e+00,  0.00000000e+00};
    //     for (int i = 0; i < size; i++)
    //         VEC1[i] -= VEC2[i];
    //     cout << sqrt(scal_prod(VEC1, VEC1)) << endl<< endl; 
    //     for (int i = 0; i < size; i++ )
    //         cout << VEC1[i] << endl;

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


double u1(int i, int j, int t)//t-номер временного слоя, ij-номера координат на сетке
{
    return sin (i * h1) 
            * sin (j * h2)
            * exp(t * tau);
}

double u2(int i, int j, int t)//t-номер временного слоя, ij-номера координат на сетке
{
    return sin (i * h1)
            * sin (j * h2)
            * exp(-t * tau);
}

double rho(int i, int j, int t)//t-номер временного слоя, ij-номера координат на сетке
{
    return (cos (i * h1) + 2) 
            * (sin (j * h2) + 2)
            * exp(t * tau);
}

double F1 (int i, int j, int t)
{
    double x1 = i * h1;
    double x2 = i * h2;
    double tt = t * tau;
    
    return 1 
        + u1(i,j,t) / rho(i,j,t) * (-sin(x1) * (sin(x2) + 2) * exp(tt))
        + u2(i,j,t) / rho(i,j,t) * ( (cos(x1) + 2) * cos(x2) * exp(tt))
        + cos(x1) * sin(x2) * exp(tt)
        + sin(x1) * cos(x2) * exp(-tt)
        ;

}

double F2 (int i, int j, int t)
{
    double x1 = i * h1;
    double x2 = i * h2;
    double tt = t * tau;


    return u1(i,j,t)
            + u1(i,j,t) * (cos(x1) * sin(x2) * exp(tt))
            + u2(i,j,t) * (sin(x1) * cos(x2) * exp(tt))
            + c / rho(i,j,t) * (-sin(x1) * (sin(x2) + 2) * exp(tt))
            - mu / rho(i,j,t) * 
                            (
                                4./3. * (-u1(i,j,t))
                                + (-u1(i,j,t))
                                +1/3. * cos(x1) * cos(x2) * exp(-tt)
                            )

            ;

}

double F3 (int i, int j, int t)
{

    double x1 = i * h1;
    double x2 = i * h2;
    double tt = t * tau;

    return -u2(i, j, t)
            + u1(i,j,t) * (cos(x1) * sin(x2) * exp(-tt))
            + u2(i,j,t) * (sin(x1) * cos(x2) * exp(-tt))
            + c / rho(i,j,t) * ((cos(x1) + 2) * cos(x2) * exp(tt))
            - mu / rho(i,j,t) * 
                            (
                                4./3. * (-u2(i,j,t))
                                + (-u2(i,j,t))
                                +1./3. * cos(x1) * cos(x2) * exp(tt)
                            )
            ;

}



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
    // if (argc != 6)
    // {
    //     cout << "USAGE: ./a.out m1 m2 n mu c\n";
    //     return -1;
    // }

    // m1 = atoi(argv[1]);
    // m2 = atoi(argv[2]);
    // n = atoi(argv[3]);
    // mu = atof(argv[4]);
    // c = atof(argv[5]);
    // Mat.resize(3*ij2k(2*m1/3, 2*m2/3) + 3);

    // tau = ((double)T) / n;
    // h1 = ((double)X1) / m1;
    // h2 = ((double)X2) / m2;

    
    // Mat.fill_matrix();
    

    Mat.resize(8);
    Mat.Mat[0].push_back(3);
    Mat.Mat[0].push_back(32);
    Mat.Mat[0].push_back(11);
    Mat.Mat[0].push_back(3);

    Mat.Ind[0].push_back(0);
    Mat.Ind[0].push_back(2);
    Mat.Ind[0].push_back(6);
    Mat.Ind[0].push_back(7);




    Mat.Mat[1].push_back(-5);
    Mat.Mat[1].push_back(-3);
    Mat.Mat[1].push_back(-1);
    Mat.Mat[1].push_back(4);

    Mat.Ind[1].push_back(0);
    Mat.Ind[1].push_back(1);
    Mat.Ind[1].push_back(6);
    Mat.Ind[1].push_back(7);

    

    Mat.Mat[2].push_back(7);
    Mat.Mat[2].push_back(1);
    Mat.Mat[2].push_back(23);
    Mat.Mat[2].push_back(29);
    Mat.Mat[2].push_back(3);
    Mat.Mat[2].push_back(4);

    Mat.Ind[2].push_back(0);
    Mat.Ind[2].push_back(2);
    Mat.Ind[2].push_back(3);
    Mat.Ind[2].push_back(4);
    Mat.Ind[2].push_back(6);
    Mat.Ind[2].push_back(7);

    Mat.Mat[3].push_back(1);
    Mat.Mat[3].push_back(12);
    Mat.Mat[3].push_back(1);
    Mat.Mat[3].push_back(32);
    Mat.Mat[3].push_back(4);
    

    Mat.Ind[3].push_back(0);
    Mat.Ind[3].push_back(1);
    Mat.Ind[3].push_back(3);
    Mat.Ind[3].push_back(6);
    Mat.Ind[3].push_back(7);



    Mat.Mat[4].push_back(-4);
    Mat.Mat[4].push_back(3);
    Mat.Mat[4].push_back(333);
    Mat.Mat[4].push_back(4);


    Mat.Ind[4].push_back(0);
    Mat.Ind[4].push_back(3);
    Mat.Ind[4].push_back(6);
    Mat.Ind[4].push_back(7);


    Mat.Mat[5].push_back(3);
    Mat.Mat[5].push_back(-3);
    Mat.Mat[5].push_back(3);
    Mat.Mat[5].push_back(4);
    Mat.Mat[5].push_back(3);
    Mat.Mat[5].push_back(23);

    Mat.Ind[5].push_back(0);
    Mat.Ind[5].push_back(1);
    Mat.Ind[5].push_back(4);
    Mat.Ind[5].push_back(5);
    Mat.Ind[5].push_back(6);
    Mat.Ind[5].push_back(7);
    


    Mat.Mat[6].push_back(4);
    Mat.Mat[6].push_back(2);
    Mat.Mat[6].push_back(132);
    Mat.Mat[6].push_back(4);
    Mat.Mat[6].push_back(1);
    Mat.Mat[6].push_back(2);


    Mat.Ind[6].push_back(0);
    Mat.Ind[6].push_back(2);
    Mat.Ind[6].push_back(3);
    Mat.Ind[6].push_back(5);
    Mat.Ind[6].push_back(6);
    Mat.Ind[6].push_back(7);



    Mat.Mat[7].push_back(2);
    Mat.Mat[7].push_back(1);
    Mat.Mat[7].push_back(3);
    Mat.Mat[7].push_back(3);
    Mat.Mat[7].push_back(3);
    Mat.Mat[7].push_back(112);
    Mat.Mat[7].push_back(2);

    Mat.Ind[7].push_back(1);
    Mat.Ind[7].push_back(2);
    Mat.Ind[7].push_back(3);
    Mat.Ind[7].push_back(4);
    Mat.Ind[7].push_back(5);
    Mat.Ind[7].push_back(6);
    Mat.Ind[7].push_back(7);


    


    Mat.b[0] = 3;
    Mat.b[1] = 2;
    Mat.b[2] = 2;
    Mat.b[3] = 21;
    Mat.b[4] = 1;
    Mat.b[5] = 1;
    Mat.b[6] = 1;
    Mat.b[7] = 1;

    vector<double> x = Mat.solveEQU();

    
    return 0;
}
