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
    int size;//размер матрицы


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

    vector<double> solveEQU()//CGS
    {
        vector<double> d(size);
        
        for (double& i: d)
            i = 1;

        vector<double> z_a(size);

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

        int MAX_ITER = 1000;
        int iter = 0;
        double alpha = 0, beta = 0;
        
        
        for (int tt = 0; tt < MAX_ITER ; tt++)
        {
            iter++;
            sparse_matr_prod_vec(p, temp);
            alpha = scal_prod(r, z_a) / scal_prod(temp, z_a);
            
            for (int i = 0; i < size; i++)
                q[i] = u[i] - alpha * temp[i];

            for (int i = 0; i < size; i++)
                d[i] = d[i] + alpha * (u[i] + q[i]);
            
            
            for (int i = 0; i < size; i++)
                temp1[i] = u[i] + q[i];
            sparse_matr_prod_vec(temp1, temp);


            for (int i = 0; i < size; i++)
                r_1[i] = r[i] - alpha * temp[i];

            beta = scal_prod(r_1, z_a) / scal_prod(r, z_a);
            for (int i = 0; i < size; i++)
                u[i] = r_1[i] + beta * q[i];
            
            for (int i = 0; i < size; i++)
                p[i] = u[i] + beta * (q[i] + beta * p[i]);

            for (int i = 0; i < size; i++)
                r[i] = r_1[i];

            if (sqrt(scal_prod(r,r)) < EPS)
                break;

        }


        sparse_matr_prod_vec (d, temp);
    
        for (int i = 0; i < size; i++)
            temp[i] -= b[i];
        double norm1 = sqrt (scal_prod(temp, temp));

        printf ("\tResidual Solve Equations: %e\n", norm1);
        printf ("\tITERATIONS: %d\n", iter);


        return d;
    }


    void fill_matrix()
    {
        int K = ij2k(m1/3*2, m2/3*2);
        vector<double> VEC1(3 * K);
        vector<double> VEC2(3 * K);
        for (int i = 0; i < K; i++)
        {
            int i1, j1;
            k2ij(i, i1, j1);
            VEC1[i] = rho(i1, j1, 0);
            VEC1[i + K] = u1(i1, j1, 0);
            VEC1[i + 2*K] = u2(i1, j1, 0);
            
        }
        int t = 0;
        int num_eq = 0; 
        
        double mu_wave = abs(VEC1[0]);
        for (int i = 0; i < K; i++)
            mu_wave = max (abs(VEC1[i]), mu_wave);
        mu_wave *= mu; 
        
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
                                        + abs( VEC1[K + ij2k(i, j)] ) / h1
                                        + abs( VEC1[2*K + ij2k(i, j)] ) / h2
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_m2


                                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                ////////////////////////////////////G_m1_m2-1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                ////////////////////////////////////G_m1-1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[K + ij2k(i, j)] - abs( VEC1[K + ij2k(i, j)] )) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                ////////////////////////////////////G_m1+1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[2*K + ij2k(i, j)] - abs( VEC1[2*K + ij2k(i, j)] )) / (2 * h2)
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
                                        + abs( VEC1[K + ij2k(i, j)] ) / h1
                                        + abs( VEC1[2*K + ij2k(i, j)] ) / h2
                );
                Ind[num_eq].push_back( ij2k(i, j) );
                //////////////////////////////////// G_m1_m2


                                
                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2 * h2)
                );
                Ind[num_eq].push_back( ij2k(i, j-1) );
                ////////////////////////////////////G_m1_m2-1



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -( abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                ////////////////////////////////////G_m1-1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[K + ij2k(i, j)] - abs( VEC1[K + ij2k(i, j)] )) / (2 * h1)
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                ////////////////////////////////////G_m1+1_m2



                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (VEC1[2*K + ij2k(i, j)] - abs( VEC1[2*K + ij2k(i, j)] )) / (2 * h2)
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

                b[num_eq] = VEC1[ij2k(i, j)] / tau  + F1(i,j,t);
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
                                        + abs( VEC1[K + ij2k(i, j)] ) / h1
                                        + abs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h1 * h1) + 2/(h2 * h2))
                );
                Ind[num_eq].push_back( K + ij2k(i, j) );
                //////////////////////////////////// V1_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j-1) );
                //////////////////////////////////// V1_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j+1) );
                //////////////////////////////////// V1_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                //////////////////////////////////// G_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                //////////////////////////////////// G_m1-1_m2


                b[num_eq] = VEC1[ij2k(i, j)] / tau 
                            + F2(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        4./3. * ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
                                        + ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
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
                                        + abs( VEC1[K + ij2k(i, j)] ) / h1
                                        + abs( VEC1[2*K + ij2k(i, j)] ) / h2
                                        + mu_wave * (4./3 * 2/(h1 * h1) + 2/(h2 * h2))
                );
                Ind[num_eq].push_back( K + ij2k(i, j) );
                //////////////////////////////////// V1_m1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j-1) );
                //////////////////////////////////// V1_m1_m2-1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -(abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i-1, j) );
                //////////////////////////////////// V1_m1-1_m2

                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-abs( VEC1[K + ij2k(i, j)] ) + VEC1[K + ij2k(i, j)]) / (2*h1)
                                        - mu_wave * 4./3. * 1./(h1*h1)
                );
                Ind[num_eq].push_back( K + ij2k(i+1, j) );
                //////////////////////////////////// V1_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        (-abs( VEC1[2*K + ij2k(i, j)] ) + VEC1[2*K + ij2k(i, j)]) / (2*h2)
                                        - mu_wave / (h2*h2)
                );
                Ind[num_eq].push_back( K + ij2k(i, j+1) );
                //////////////////////////////////// V1_m1_m2+1


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        c
                );
                Ind[num_eq].push_back( ij2k(i+1, j) );
                //////////////////////////////////// G_m1+1_m2


                ////////////////////////////////////
                Mat[num_eq].push_back(
                                        -c
                );
                Ind[num_eq].push_back( ij2k(i-1, j) );
                //////////////////////////////////// G_m1-1_m2


                b[num_eq] = VEC1[ij2k(i, j)] / tau 
                            + F2(i,j,t)
                            - (mu_wave - mu * exp( -VEC1[ij2k(i, j)] ) )
                                    * (
                                        4./3. * ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
                                        + ( VEC1[K + ij2k(i-1, j)] - 2 * VEC1[K + ij2k(i, j)] + VEC1[K + ij2k(i+1, j)]) / (h1*h1)
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
        ///////////////////////////////////////////////////////////////////////////////////////////////// 3-е уравнение
        for (int i = 1; i < m1/3; i++)
        {
            for (int j = 1; j < m2; j++)
            {
                Mat[num_eq].clear();
                Ind[num_eq].clear();
                

                num_eq++;


            }
        }

        for (int i = 0; i < m1/3; i++)
        {
            for (int j = 1; j < m2/3; j++)
            {
                i += m1/3;
                j += m2/3;


                num_eq++;
                i -= m1/3;
                j -= m2/3;
            }
        }
        ///////////////////////////////////////////////////////////////////////////////////////////////////////////////





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
    return sin (2 * M_PI * i * h1) * 
            sin (2 * M_PI * j * h2)
            * exp(t * tau);
}

double u2(int i, int j, int t)//t-номер временного слоя, ij-номера координат на сетке
{
    return sin (2 * M_PI * i * h1) * 
            sin (2 * M_PI * j * h2)
            * exp(-t * tau);
}

double rho(int i, int j, int t)//t-номер временного слоя, ij-номера координат на сетке
{
    return (cos (2 * M_PI * i * h1) + 1.5) *
            (sin (2 * M_PI * j * h2) + 1.5)
            * exp(t * tau);
}

double F1 (int i, int j, int t)
{
    return 1 
    + u1(i, j, t) * 1. / rho(i, j, t) * exp(t * tau) * (sin(2 * M_PI * j * h2) + 1.5) * (-2 * M_PI * sin(2 * M_PI * h1 * i))
    + u2(i, j, t) * 1. / rho(i, j, t) * exp(t * tau) * (cos(2 * M_PI * i * h1) + 1.5) * (2 * M_PI * cos(2 * M_PI * h2 * j))
    + sin (2 * M_PI * h2 * j) * cos(2 * M_PI * h1 * i) * 2 * M_PI * exp(t * tau)
    + sin (2 * M_PI * h1 * i) * cos(2 * M_PI * h2 * j) * 2 * M_PI * exp(-t * tau);
}

double F2 (int i, int j, int t)
{
    return u1(i, j, t)
            + u1(i,j,t) * exp(t * tau) * sin(2 * M_PI * h2 * j) * cos (2 * M_PI * h1 * i) * 2 * M_PI
            + u2(i,j,t) * exp(t * tau) * sin(2 * M_PI * h1 * i) * cos (2 * M_PI * h2 * j) * 2 * M_PI
            - mu/rho(i,j,t) * (
                            4./3. * (- 4 * M_PI * M_PI * u1(i,j,t)) 
                            + (-4 * M_PI * M_PI * u1(i,j,t))
                            + 1/3. * 4 * M_PI * M_PI * cos(2 * M_PI * h1 * i) * cos(2 * M_PI * h2 * j) * exp(-t * tau)
                        )
            + c * exp(t) * (sin(2*M_PI * h2 * j) + 1.5) * (-sin(2*M_PI * h1 * i) * 2 * M_PI) * 1./rho(i,j,t);
            ;

}

double F3 (int i, int j, int t)
{
    return u2(i, j, t)
            + u1(i,j,t) * exp(-t * tau) * sin(2 * M_PI * h2 * j) * cos (2 * M_PI * h1 * i) * 2 * M_PI
            + u2(i,j,t) * exp(-t * tau) * sin(2 * M_PI * h1 * i) * cos (2 * M_PI * h2 * j) * 2 * M_PI
            - mu/rho(i,j,t) * (
                            4./3. * (- 4 * M_PI * M_PI * u2(i,j,t)) 
                            + (-4 * M_PI * M_PI * u2(i,j,t))
                            + 1/3. * 4 * M_PI * M_PI * cos(2 * M_PI * h1 * i) * cos(2 * M_PI * h2 * j) * exp(t * tau)
                        )
            + c * exp(t) * (cos(2*M_PI * h1 * i) + 1.5) * (cos(2*M_PI * h2 * j) * 2 * M_PI) * 1./rho(i,j,t)
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
    if (k <= (m2 + 1) * (m1 / 3 + 1))
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
    if (argc != 6)
    {
        cout << "wrong argc\n";
        return -1;
    }

    m1 = atoi(argv[1]);
    m2 = atoi(argv[2]);
    n = atoi(argv[3]);
    mu = atof(argv[4]);
    c = atof(argv[5]);
    Mat.resize(3*ij2k(2*m1/3, 2*m2/3));

    
    return 0;
}
