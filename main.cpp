#include<stdio.h>
#include<iostream>
#include<vector>
using namespace std;

int m1, m2;
int T=1, X1=1, X2=1;


class Matrix
{
public:
    vector<vector<double>> Mat;
    vector<vector<int>> Ind;
    int size;

    Matrix (int n): size(n)
    {
        Mat.resize(n);
        Ind.resize(n);
    }

};

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

int main()
{
    Matrix a(20);
    m1=m2=9;
    int i,j,k=47;
    k2ij(k,i,j);
    cout << i << " " << j << endl;
    printf("hhahah\n");
    return 0;
}
