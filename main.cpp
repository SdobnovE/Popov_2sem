#include<stdio.h>
#include<iostream>
#include<vector>
#include<math.h>
#include <stdlib.h>
using namespace std;

const double EPS = 1e-16;

int m1, m2, n;
double mu;
int T=1, X1=1, X2=1;
double tau, h1, h2;


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
        double alpha = 0, beta = 0;
        
        
        for (int tt = 0; tt < MAX_ITER ; tt++)
        {
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


        return d;
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
    return 1 + u1(i, j, t) * 1. / rho(i, j, t) * exp(t * tau) * (sin(2 * M_PI * j * h2) + 1.5) * (-2 * M_PI * sin(2 * M_PI * h1 * i))
    + u2(i, j, t) * 1. / rho(i, j, t) * exp(t * tau) * (cos(2 * M_PI * i * h1) + 1.5) * (2 * M_PI * cos(2 * M_PI * h2 * j))
    + sin (2 * M_PI * h2 * j) * cos(2 * M_PI * h1 * i) * 2 * M_PI * exp(t * tau)
    + sin (2 * M_PI * h1 * i) * cos(2 * M_PI * h2 * j) * 2 * M_PI * exp(-t * tau);
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
    if (argc != 5 )
    {
        cout << "wrong argc\n";
        return -1;
    }
    // m1 = atoi(argv[1]);
    // m2 = atoi(argv[2]);
    // n = atoi(argv[3]);
    // mu = atof(argv[4]);
    // Mat.resize(ij2k(2*m1/3, 2*m2/3));

    Mat.resize(5);
    Mat.Mat[0].push_back(11);
    Mat.Mat[0].push_back(3);

    Mat.Mat[1].push_back(1);
    
    Mat.Mat[2].push_back(1);
    
    Mat.Mat[3].push_back(112);
    Mat.Mat[3].push_back(2);

    Mat.Mat[4].push_back(2);
    Mat.Mat[4].push_back(1);




    Mat.Ind[0].push_back(0);
    Mat.Ind[0].push_back(2);

    Mat.Ind[1].push_back(1);
    
    Mat.Ind[2].push_back(2);
    
    Mat.Ind[3].push_back(0);
    Mat.Ind[3].push_back(3);
    
    Mat.Ind[4].push_back(2);
    Mat.Ind[4].push_back(4);

    Mat.b[0] = 13;
    Mat.b[1] = 14;
    Mat.b[2] = 13;
    Mat.b[3] = 14;
    Mat.b[4] = 13;

    vector<double> x = Mat.solveEQU();
    
    for (auto i : x)
        cout << i << endl;
    return 0;
}
