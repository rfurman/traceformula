#define NDEBUG
#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <iomanip>
#include <Eigen/MPRealSupport>

using namespace std;
using Eigen::Matrix;
using Eigen::Dynamic;

typedef mpfr::mpreal scalar;
//typedef complex<mpfr::mpreal> scalar;

int main(int argc, char **argv) {
    cout << "Testing eigenvalues" << endl;
    int MM, prec;
    cin >> MM >> prec;
    mpfr::mpreal::set_default_prec(prec);

    cout << "Random matrices" << endl;
    Matrix<scalar,Dynamic,Dynamic> mat1(MM,MM),mat2(MM,MM);
    for(int i=0; i<MM; i++) for(int j=0; j<=i; j++) {
        mat1(i,j) = mat1(j,i) = rand()/65536./65536.;
        mat2(i,j) = mat2(j,i) = rand()/65536./65536.;
    }
    cout << "Product" << endl;
    //Matrix<scalar,Dynamic,Dynamic> mat=mat2.inverse()*mat1;

    if(argv[1][0]=='D') {
        cout << "Diagonalisation" << endl;
        //Eigen::ComplexEigenSolver<Matrix<scalar,Dynamic,Dynamic> > ces(mat);
        //cout << setprecision(0) << ces.eigenvalues() << endl;
    }

    if(argv[1][0]=='G') {
        cout << "General Diagonalisation" << endl;
        //Eigen::GeneralizedEigenSolver<Matrix<scalar,Dynamic,Dynamic> > ces(mat1, mat2);
        //cout << setprecision(0) << ces.eigenvalues() << endl;
    }

    if(argv[1][0]=='S') {
        cout << "General Diagonalisation" << endl;
        Eigen::GeneralizedSelfAdjointEigenSolver<Matrix<scalar,Dynamic,Dynamic> > ces(mat1, mat2);
        //cout << setprecision(0) << ces.eigenvalues() << endl;
    }


    return 0;
}
