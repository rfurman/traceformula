// build with something like
//
// g++ arbtest.cc -fpermissive -I. -O2 -lflint -o arbtest -I/home/bober/include/flint

#include <iostream>
#include <cstdlib>
#include <ctime>
#include <complex>
#include <arb/fmpcb_mat.h>

#include "Eigen/Dense"

using std::cout;
using std::endl;
using std::flush;
using std::complex;

using Eigen::Matrix;
using Eigen::ComplexEigenSolver;
using Eigen::Dynamic;

typedef Matrix<complex<double>, Dynamic, Dynamic> complex_matrix_t;

double seconds_since_last_call() {
    static clock_t lasttime = 0;
    clock_t thistime = clock();
    double seconds = (thistime - lasttime)/(double)CLOCKS_PER_SEC;
    lasttime = thistime;
    return seconds;
}

void fmpcb_set_complex(fmpcb_t z, complex<double> w, double error = 1e-15) {
    fmpr_set_d(fmprb_midref(fmpcb_realref(z)), w.real());
    fmpr_set_d(fmprb_midref(fmpcb_imagref(z)), w.imag());
    fmpr_set_d(fmprb_radref(fmpcb_realref(z)), error);
    fmpr_set_d(fmprb_radref(fmpcb_imagref(z)), error);
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


    cout << "Converting to arb... " << flush;
    int prec = 750;

    fmpcb_mat_t Q, Qinv, AA, L;
    fmpcb_mat_init(Q, N, N);
    fmpcb_mat_init(Qinv, N, N);
    fmpcb_mat_init(AA, N, N);
    fmpcb_mat_init(L, N, N);

    for(int j = 0; j < N; j++) {
        for(int k = 0; k < N; k++) {
            fmpcb_set_complex(fmpcb_mat_entry(AA, j,k), A(j,k), 0);
            fmpcb_set_complex(fmpcb_mat_entry(Q, j,k), ces.eigenvectors()(j,k), 0);
        }
    }
    cout << seconds_since_last_call() << endl;

    cout << "Doing some multiplications." << endl;
    cout << "inverting... " << flush;
    if(!fmpcb_mat_inv(Qinv, Q, prec)) {
        cout << "Matrix could not be inverted" << endl;
        return -1;
    }
    cout << seconds_since_last_call() << endl;

    fmpcb_printd(fmpcb_mat_entry(Qinv, 0, 0), 10);
    cout << endl;
    
    cout << "multiplying..." << flush;
    fmpcb_mat_mul(L, Qinv, AA, prec);
    cout << seconds_since_last_call() << endl << "multiplying..." << flush;
    fmpcb_mat_mul(L, L, Q, prec);
    cout << seconds_since_last_call() << endl;

    fmpcb_printd(fmpcb_mat_entry(L,0,0), 10);
    cout << " " << LL(0,0) << endl;
    fmpcb_printd(fmpcb_mat_entry(L,0,1), 10);
    cout << " " << LL(0,1) << endl;
    fmpcb_printd(fmpcb_mat_entry(L,0,2), 10);
    cout << " " << LL(0,2) << endl;

    return 0;
}
