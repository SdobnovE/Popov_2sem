#include<stdio.h>
#include<iostream>
#include<vector>
#include <stdlib.h>
using namespace std;

int m1, m2, n;
double mu;
int T=1, X1=1, X2=1;

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
    vector<vector<double>> Mat;
    vector<vector<int>> Ind;
    int size;


    Matrix();

    Matrix (int n): size(n)
    {
        Mat.resize(n);
        Ind.resize(n);
    }
    
    void resize (int n)
    {
        size = n;
        Mat.resize(n);
        Ind.resize(n);
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
    m1 = atoi(argv[1]);
    m2 = atoi(argv[2]);
    n = atoi(argv[3]);
    mu = atof(argv[4]);
    Mat.resize(ij2k(2*m1/3, 2*m2/3));

    
    
    return 0;
}
