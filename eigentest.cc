// build with something like
// g++ eigentest.cc -I. -O2 -lmpfr -o eigentest

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <complex>

#include "Eigen/Dense"
#include "Eigen/MPRealSupport"

using std::cout;
using std::endl;
using std::flush;
using std::complex;

using Eigen::Matrix;
using Eigen::ComplexEigenSolver;
using Eigen::Dynamic;

typedef Matrix<complex<double>, Dynamic, Dynamic> complex_matrix_t;
typedef mpfr::mpreal RR;
typedef complex<RR> CC;
typedef Matrix<CC,Dynamic,Dynamic> CCMatrix;
typedef Matrix<RR,Dynamic,Dynamic> RRMatrix;


double seconds_since_last_call() {
    static clock_t lasttime = 0;
    clock_t thistime = clock();
    double seconds = (thistime - lasttime)/(double)CLOCKS_PER_SEC;
    lasttime = thistime;
    return seconds;
}

int main(int argc, char ** argv) {
    if(argc > 1)
        std::srand(std::atoi(argv[1]));
    else
        std::srand(std::time(NULL));

    int N = 250;

    cout << "Building matrix." << endl;

    complex_matrix_t A(N,N);
    for(int j = 0; j < N; j++) {
        for(int k = 0; k < N; k++) {
            A(j,k) = std::rand();
        }
    }

    cout << "Computing eigenvalues... " << flush;
    ComplexEigenSolver<complex_matrix_t> ces;
    ces.compute(A, true);
    //cout << ces.eigenvalues() << endl;
    //cout << ces.eigenvectors() << endl;

    complex_matrix_t LL = (ces.eigenvectors().inverse()) * A * ces.eigenvectors();
    cout << seconds_since_last_call() << endl;


    cout << "Converting to high precision... " << flush;
    int prec = 750;
    mpfr::mpreal::set_default_prec(prec);

    CCMatrix Q(N,N), Qinv(N,N), AA(N,N), L(N,N);

    for(int j = 0; j < N; j++) {
        for(int k = 0; k < N; k++) {
            AA(j,k) = A(j,k);
            Q(j,k) = ces.eigenvectors()(j,k);
        }
    }
    cout << seconds_since_last_call() << endl;

    cout << "Doing some multiplications." << endl;
    cout << "inverting... " << flush;
    Qinv = Q.inverse();
    cout << seconds_since_last_call() << endl;

    cout << Qinv(0,0) << endl;
    
    cout << "multiplying..." << flush;
    L = Qinv * AA;
    cout << seconds_since_last_call() << endl << "multiplying..." << flush;
    L = L * Q;
    cout << seconds_since_last_call() << endl;

    cout << L(0,0) << " " << LL(0,0) << endl;
    cout << L(0,1) << " " << LL(0,1) << endl;
    cout << L(0,2) << " " << LL(0,2) << endl;

    return 0;
}
