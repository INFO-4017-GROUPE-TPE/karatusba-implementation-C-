/*============================================================================
// Name        : KaratsubaMultiplication
// Author      : 
============================================================================

	| NOUTCHEU LIBERT | 19M2310
	|
	|
	|
	|
	|

============================================================================
// Version     : 1.0
// Description : Programme en C++  Matrice Multiplications avec la methode de  Strassen 
//============================================================================
*/


#include <bits/stdc++.h>
#include <cmath>

using namespace std;
#define table vector<int>
#define matrice vector<table>
//! addition and subtraction
void add(matrice &A, matrice &B, matrice &C, int size);
void ConvertToSquareMat(matrice &A, matrice &B, int r1, int c1, int r2, int c2);
void reduce_matrix(matrice &A, matrice &B, matrice &C, int size);

void Strassen_algorithm(matrice &A, matrice &B, matrice &C, int size);

/* recherche la prochaine puissance de 2 */
int nextPowerOf2(int k)
{
    return pow(2, int(ceil(log2(k))));
}

// affice le resultat de la  matrix
void display(matrice C, int m, int n)
{
    for (int i = 0; i < m; i++)
    {
        cout << "|"
             << " ";
        for (int j = 0; j < n; j++)
        {
            cout << C[i][j] << " ";
        }
        cout << "|" << endl;
    }
}

// AJOUTER LES DONNES DANS UN MATRICE CARRE
void insert(matrice& x, int n)
{
	double val;
	for(int i=0;i<n; i++)
	{
		for(int j=0;j<n;j++)
		{
			cin>>val;
			x[i][j]=val;
		}
	}
}




int main()
{   
    int matrix_size = 0;
    cout<<"\n define the order of a matrix \n";
    cin>>matrix_size;
    table z(matrix_size);

    matrice a(matrix_size, z) ;
    matrice b(matrix_size, z) ;

    cout<<"\n a: \n";
	insert(a,matrix_size);
	//insert values in the matrix a
	cout<<"\n b: \n";
	insert(b, matrix_size);
    cout<<"\n define the order of a matrix \n";
    ConvertToSquareMat(a, b, matrix_size, matrix_size, matrix_size, matrix_size); // A[][],B[][],R1,C1,R2,C2
    return 0;
}


//! addition and subtraction
void add(matrice &A, matrice &B, matrice &C, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] + B[i][j];
        }
    }
}



void reduce_matrix(matrice &A, matrice &B, matrice &C, int size)
{
    for (int i = 0; i < size; i++)
    {
        for (int j = 0; j < size; j++)
        {
            C[i][j] = A[i][j] - B[i][j];
        }
    }
}


/*pour transformer la matrix en matrice carre */

void ConvertToSquareMat(matrice &A, matrice &B, int r1, int c1, int r2, int c2)
{
    int maxSize = max({r1, c1, r2, c2});
    int size = nextPowerOf2(maxSize);
 
    table z(size);
    matrice Aa(size, z), Bb(size, z), Cc(size, z);
 
    for (unsigned int i = 0; i < r1; i++)
    {
        for (unsigned int j = 0; j < c1; j++)
        {
            Aa[i][j] = A[i][j];
        }
    }
    for (unsigned int i = 0; i < r2; i++)
    {
        for (unsigned int j = 0; j < c2; j++)
        {
            Bb[i][j] = B[i][j];
        }
    }

    Strassen_algorithm(Aa, Bb, Cc, size);
    table resutl_matrix(c2);
    matrice C(r1, resutl_matrix);
    for (unsigned int i = 0; i < r1; i++)
    {
        for (unsigned int j = 0; j < c2; j++)
        {
            C[i][j] = Cc[i][j];
        }
    }
    display(C, r1, c1);
}



//!-----------------------------
void Strassen_algorithm(matrice &A, matrice &B, matrice &C, int size)
{
    if (size == 1)
    {
        C[0][0] = A[0][0] * B[0][0];
        return;
    }
    else
    {
        int newSize = size / 2;
        table z(newSize);
        matrice a(newSize, z), b(newSize, z), c(newSize, z), d(newSize, z),
            e(newSize, z), f(newSize, z), g(newSize, z), h(newSize, z),
            c11(newSize, z), c12(newSize, z), c21(newSize, z), c22(newSize, z),
            p1(newSize, z), p2(newSize, z), p3(newSize, z), p4(newSize, z),
            p5(newSize, z), p6(newSize, z), p7(newSize, z), fResult(newSize, z),
            sResult(newSize, z);
        int i, j;
 
        //! diviser la matrice ne parts equitable
        for (i = 0; i < newSize; i++)
        {
            for (j = 0; j < newSize; j++)
            {
                a[i][j] = A[i][j];
                b[i][j] = A[i][j + newSize];
                c[i][j] = A[i + newSize][j];
                d[i][j] = A[i + newSize][j + newSize];
 
                e[i][j] = B[i][j];
                f[i][j] = B[i][j + newSize];
                g[i][j] = B[i + newSize][j];
                h[i][j] = B[i + newSize][j + newSize];
            }
        }
        /*
             A         B           C
           [a b]   * [e f]   =  [c11 c12]
                [g h]      [c21 c22]
           p1,p2,p3,p4=AHED for this: A:Row(+) and B:Column(-)
           p5=Diagonal :both +ve
           p6=Last CR  :A:Row(-) B:Column(+)
           p7=First CR :A:Row(-) B:Column(+)
        */
        //! calcule de tous les formule de  strassen 
        //*p1=a*(f-h)
        reduce_matrix(f, h, sResult, newSize);
        Strassen_algorithm(a, sResult, p1, newSize);
 
        //*p2=h*(a+b)
        add(a, b, fResult, newSize);
        Strassen_algorithm(fResult, h, p2, newSize);
 
        //*p3=e*(c+d)
        add(c, d, fResult, newSize);
        Strassen_algorithm(fResult, e, p3, newSize);
 
        //*p4=d*(g-e)
        reduce_matrix(g, e, sResult, newSize);
        Strassen_algorithm(d, sResult, p4, newSize);
 
        //*p5=(a+d)*(e+h)
        add(a, d, fResult, newSize);
        add(e, h, sResult, newSize);
        Strassen_algorithm(fResult, sResult, p5, newSize);
 
        //*p6=(b-d)*(g+h)
        reduce_matrix(b, d, fResult, newSize);
        add(g, h, sResult, newSize);
        Strassen_algorithm(fResult, sResult, p6, newSize);
 
        //*p7=(a-c)*(e+f)
        reduce_matrix(a, c, fResult, newSize);
        add(e, f, sResult, newSize);
        Strassen_algorithm(fResult, sResult, p7, newSize);
 
        /* calucle de tous les element de  C avec p1,p2,p3, p4
        c11=p4+p5+p6-p2
        c12=p1+p2
        c21=p3+p4
        c22=p1-p3+p5-p7
        */
        add(p1, p2, c12, newSize); //!
        add(p3, p4, c21, newSize); //!
 
        add(p4, p5, fResult, newSize);
        add(fResult, p6, sResult, newSize);
        reduce_matrix(sResult, p2, c11, newSize); //!
 
        reduce_matrix(p1, p3, fResult, newSize);
        add(fResult, p5, sResult, newSize);
        reduce_matrix(sResult, p7, c22, newSize); //!
 
        // regrouper le resultat obtenu en une seul matrix:
        for (i = 0; i < newSize; i++)
        {
            for (j = 0; j < newSize; j++)
            {
                C[i][j] = c11[i][j];
                C[i][j + newSize] = c12[i][j];
                C[i + newSize][j] = c21[i][j];
                C[i + newSize][j + newSize] = c22[i][j];
            }
        }
    }
}
