#include <Eigen/Dense>
#include <complex>
#include <iostream>
#include <iomanip>
#include <Eigen/MPRealSupport>

using namespace std;
using Eigen::Matrix;
using Eigen::ComplexEigenSolver;
using Eigen::Dynamic;

int main(void) {
    cout << "Testing eigenvalues" << endl;
    int MM, prec;
    cin >> MM >> prec;
    mpfr::mpreal::set_default_prec(prec);

    Matrix<complex<mpfr::mpreal>,Dynamic,Dynamic> mat(MM,MM);
    for(int i=0; i<MM; i++) for(int j=0; j<MM; j++)
        mat(i,j) = rand()/65536./65536.;
    ComplexEigenSolver<Matrix<complex<mpfr::mpreal>,Dynamic,Dynamic> > ces(mat);
    cout << setprecision(0) << ces.eigenvalues() << endl;
    return 0;
}
