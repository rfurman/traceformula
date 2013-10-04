#include <Eigen/Dense>
#include <complex>
#include <fstream>
#include <iomanip>
#include <pari/pari.h>
#include <algorithm>
#include <iostream>
#include <cassert>
#include <vector>
#include <Eigen/MPRealSupport>

//#include <boost/multiprecision/gmp.hpp> 
//#include <boost/multiprecision/cpp_int.hpp> 
//#include <flint/fmprb_mat.h>
using namespace std;
//using namespace mpfr;
using Eigen::Matrix;
using Eigen::ComplexEigenSolver;
using Eigen::Dynamic;
using Eigen::ColMajor;
//using namespace boost::multiprecision;
typedef long long i64;
//typedef mpz_int i64;

//typedef double Real;

typedef mpfr::mpreal RR;
typedef complex<RR> CC;

typedef Matrix<CC,Dynamic,Dynamic> CCMatrix;
typedef Matrix<RR,Dynamic,Dynamic> RRMatrix;
typedef Matrix<long long,Dynamic,Dynamic,ColMajor> ZZMatrix;

#include "number_theoretic.h"


double to_double(long long x) { return x; }
//double to_double(const mpz_int& x) { return x.convert_to<double>(); }

#include "characters.h"
#include "trace_formula.h"

//#define pow boost::multiprecision::pow

// Divide through by coefficient of 1
CCMatrix adjust(CCMatrix mat) {
    for(int r=0; r<mat.rows(); r++) {
        for(int c=2; c<mat.cols(); c++) {
            mat(r,c) /= mat(r,1);
        }
        mat(r,1)=1;
    }
    return mat;
}

// Compute a(m)*a(n) via multiplicativity relations with character chi

template<class T, class U> CC prod2(const T& a, const U& chi, int m, int n, int k) {
    int g = __gcd(m,n);
    CC ret;
    for(int d=1; d<=g; d++) if(g%d==0) {
        int dk=1; // d^(k-1)
        for(int i=0; i<k-1; i++) dk*=d;
        ret += a(m*n/d/d)*RR(dk)*chi(d%chi.cols());
    }
    return ret;
} 

template<class T> bool close_enough(const T& A, const T& B) {
    return abs(A-B) < 0.00001 * (abs(A)+abs(B)+1);
}

int main(void) {
        pari_init(8000000, 600000);

        cout << "\rInitalizing primes" << flush;
        for(int i=1; i<41541; i++) primes_list.push_back(itos(prime(i)));

        cout << "\rComputing mobius mu" << flush;
        for(int i=1; i<500000; i++) mobius[i]=mobius2[i]=1;
        int p;
        for(int i=0;; i++) {
            int p = primes_list[i];
            if(p>500000) break;
            for(int i=p; i<500000; i+=p) mobius[i]*=-1;
            for(int i=p; i<500000; i+=p) mobius2[i]*=-2;
        }
        for(int i=0;; i++) {
            int p = primes_list[i];
            if(p*p>500000) break;
            for(int i=p*p; i<500000; i+=p*p) mobius[i]*=0;
            for(int i=p*p; i<500000; i+=p*p) mobius2[i]/=-2;
        }
        for(int i=0;; i++) {
            int p = primes_list[i];
            if(p*p*p>500000) break;
            for(int i=p*p*p; i<500000; i+=p*p*p) mobius2[i]*=0;
        }

        //int M=500000, N=997, k=4;
        //int M=2, N=997, k=4;
        int M=1000, N=15, k=4, prec=1000;

        cin >> M >> N >> k >> prec;
        mpfr::mpreal::set_default_prec(prec);

        int phiN=phi(N);

        if(theoretical_TrT_bound(M,N,k)<1e18) {
            // Compute traces
            RRMatrix vals = allTrThat12new(M,N,k);
            CCMatrix fourier = matrix_of_characters_by_parity(N,k);
            cout << "\rFourier transform              " << flush;
            CCMatrix vals2 = fourier * vals;
            cout << "\rDivision                       " << flush;
            vals2 /= CC(12*phiN);

            for(int chi=0; chi<fourier.rows(); chi++) {
                int dim = round(vals2(chi,1).real()).toLong();
                assert(abs(vals2(chi,1)-(RR)dim)<0.0001);
                int number_good=0, number_bad=0;
                if(dim==0) {
                    cout << "\rComputing dimension " << dim << " eigenbasis with character #" << chi+1 << endl;
                    continue;
                } else {
                    cout << "\rComputing dimension " << dim << " eigenbasis with character #" << chi+1 << endl;
                }

                vector<int> rel_primes;
                for(int i=1; rel_primes.size()<dim; i++) if(__gcd(i,2)==1) rel_primes.push_back(i);

                CCMatrix c(dim,dim), Tp(dim,dim), bb(dim, M/rel_primes.back());
                cout << "\rLoading small data                                " << flush;
                for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) c(i,j)=prod2(vals2.row(chi),fourier.row(chi),rel_primes[i],rel_primes[j],k);
                for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) Tp(i,j)=prod2(vals2.row(chi),fourier.row(chi),rel_primes[i],rel_primes[j]*2,k);
                cout << "\rLoading big data                  " << flush;
                for(int i=0; i<dim; i++) for(int j=0; j<M/rel_primes.back(); j++) bb(i,j)=prod2(vals2.row(chi),fourier.row(chi),rel_primes[i],j,k);
                cout << "\rBasic lin alg                  " << flush;
                CCMatrix v = (c.inverse() * Tp);
                cout << "\rSolve for eigenvectors           " << flush;
                ComplexEigenSolver<CCMatrix > ces((c.inverse() * Tp));
                cout << "\rCompute result                        " << flush;
                CCMatrix results = adjust(ces.eigenvectors().transpose() * bb);
                //cout << results << endl;
                //cout << ces.eigenvalues() << endl;

                cout << "\rVerifying multiplicativity             " << endl;
                for(int i=0; i<results.rows(); i++) {
                    bool good = true;
                    for(int j=0; j<results.rows(); j++) if(i!=j)
                        if(close_enough(ces.eigenvalues()(i),ces.eigenvalues()(j)))
                            good = false;
                    if(!good) {
                        cout << "\rSkipping high multiplicity eigenvalue " << ces.eigenvalues()(i) << endl;
                        number_bad++;
                        continue;
                    }
                    for(int j=1; j<results.cols(); j++)
                      for(int l=1; j*l<results.cols(); l++) {
                        CC A = results(i,j)*results(i,l);
                        CC B = prod2(results.row(i),fourier.row(chi),j,l,k);
                        if(!close_enough(A,B)) {
                            if(good) {
                                //cout << chi << " " << i << " " << j << " " << k << " " << A << " " << B << " " << results(i,j) << " " << results(i,k) << " " << results(i,j*k) << endl;
                                //cout << fourier.row(chi) << endl;
                                //assert(false);
                            }
                            good=false;
                        }
                    }
                    if(!good) {
                        cout << "\rFailed multiplicativity check        " << flush;
                        number_bad++;
                    } else {
                        cout << "\rPassed multiplicativity check         " << flush;
                        number_good++;
                    }
                    char filename[100];
                    sprintf(filename, "newforms/level%d_character%d_weight%d_number%d", N, chi, k, i);
                    ofstream out(filename);
                    for(int j=1; j<results.cols(); j++)
                        out << setprecision(prec*0.28) << results(i,j).real() << " + " << results(i,j).imag() << " i" << endl;
                    out.close();
                }

                cout << "\r" << M/rel_primes.back() << " coefficients for " << number_good << " newforms (" << number_bad << " failed to be computed)" << endl;
            }

        } else {
            cout << "Calculations would overflow long long, switching to gmp" << endl;
            //allTrThat12<mpz_int>(M,N,k);
            assert(false);
        }
}
