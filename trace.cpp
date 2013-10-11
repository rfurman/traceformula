#define NDEBUG
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
typedef long i64;
//typedef mpz_int i64;

//typedef double Real;

typedef long ZZ;
typedef mpfr::mpreal RR;
typedef complex<RR> CC;


// I thought this would fix RRMatrix*ZZMatrix, but it is not enough
namespace Eigen {
namespace internal {
template<> struct scalar_product_traits<RR,ZZ>
{
  enum {
    // Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef RR ReturnType;
};

template<> struct scalar_product_traits<CC,ZZ>
{
  enum {
    // Cost = 2*NumTraits<T>::MulCost,
    Defined = 1
  };
  typedef CC ReturnType;
};

}
}

CC operator*(CC A, const ZZ& B) {
    A *= B;
    return A;
}

typedef Matrix<CC,Dynamic,Dynamic> CCMatrix;
typedef Matrix<RR,Dynamic,Dynamic> RRMatrix;
typedef Matrix<long,Dynamic,Dynamic,ColMajor> ZZMatrix;

template<class T> CCMatrix explicit_product(const T& fourier_matrix, const ZZMatrix& mat) {
    CCMatrix ret = CCMatrix::Zero(fourier_matrix.rows(), mat.cols());
    assert(fourier_matrix.cols() == mat.rows());
    for(int i=0; i<fourier_matrix.rows(); i++)
        for(int j=0; j<mat.rows(); j++)
            for(int k=0; k<mat.cols(); k++)
                ret(i,k) += fourier_matrix(i,j)*mat(j,k);
    return ret;
}



#include "number_theoretic.h"


double to_double(long x) { return x; }
//double to_double(const mpz_int& x) { return x.convert_to<double>(); }

#include "characters.h"
#include "trace_formula.h"

//#define pow boost::multiprecision::pow

// Divide through by coefficient of 1
void adjust(CCMatrix& mat) {
    for(int r=0; r<mat.rows(); r++) {
        for(int c=2; c<mat.cols(); c++) {
            mat(r,c) /= mat(r,1);
        }
        mat(r,1)=1;
    }
}

// Compute a(m)*a(n) via multiplicativity relations with character chi

template<class T, class U> CC prod2(const T& a, const U& chi, int m, int n, int k) {
    int g = __gcd(m,n);
    CC ret;
    for(int d=1; d<=g; d++) if(g%d==0) {
        i64 dk=1; // d^(k-1)
        for(int i=0; i<k-1; i++) dk*=d;
        ret += a(m/d*n/d)*RR(dk)*chi(d%chi.cols());
    }
    return ret;
} 

template<class T, class U> CC prod3(const T& a, const U& chi, int m, int n, int p, int k) {
    int g = __gcd(m,n);
    CC ret;
    for(int d=1; d<=g; d++) if(g%d==0) {
        i64 dk=1; // d^(k-1)
        for(int i=0; i<k-1; i++) dk*=d;
        ret += prod2(a,chi,m/d*n/d,p,k)*RR(dk)*chi(d%chi.cols());
    }
    return ret;
} 



template<class T> bool close_enough(const T& A, const T& B) {
    return abs(A-B) < 1e-100 * (abs(A)+abs(B)+1);
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
            ZZMatrix vals = allTrThat12new(M,N,k);
            //for(int a=0; a<vals.rows(); a++)
                //for(int b=0; b<vals.cols(); b++)
                    //cout << a<< " " << b << " : " << vals(a,b) << endl;
            CCMatrix fourier = matrix_of_characters_by_parity(N,k);

            for(int chi=0; chi<fourier.rows(); chi++) {
                cout << "\rFourier transform              " << flush;
                CCMatrix vals2 = explicit_product(fourier.row(chi), vals);
                cout << "\rDivision                       " << flush;
                vals2 /= CC(12*phiN);

                int dim = round(vals2(0,1).real()).toLong();
                if(abs(vals2(0,1)-(RR)dim)>0.0001) abort();
                int number_good=0, number_bad=0;
                if(dim==0) {
                    cout << "\rComputing dimension " << dim << " eigenbasis with character #" << chi+1 << endl;
                    continue;
                } else {
                    cout << "\rComputing dimension " << dim << " eigenbasis with character #" << chi+1 << endl;
                }

                cout << "\rPicking c matrix                                " << flush;
                vector<int> rel_primes;
                CCMatrix c(dim,dim);

                for(int i=1; rel_primes.size()<dim; i++) if(__gcd(i,N)==1) {
                    rel_primes.push_back(i);
                    int j=rel_primes.size()-1;
                    for(int i=0; i<rel_primes.size(); i++) c(i,j)=c(j,i)=prod2(vals2.row(0),fourier.row(chi),rel_primes[i],rel_primes[j],k);
                }
 
                if(abs(c.determinant())<0.000001) {
                    // If the candidate matrix is not invertible, build it up instead more conservatively
                    //    keeping full rank at each step.
                    rel_primes.clear();
                    for(int i=1; rel_primes.size()<dim; i++) if(__gcd(i,N)==1) {
                        rel_primes.push_back(i);
                        int j=rel_primes.size()-1;
                        for(int i=0; i<rel_primes.size(); i++) c(i,j)=c(j,i)=prod2(vals2.row(0),fourier.row(chi),rel_primes[i],rel_primes[j],k);
                        if(abs(c.topLeftCorner(j+1,j+1).determinant())<0.000001) rel_primes.pop_back();
                    }
                    cout << "\rHad to skip " << rel_primes.back()-dim << " elements to get a rank " << dim << " c-matrix" << endl;
                }
                

                // Number of sqrt(chi) to generate
                int Nchi = max(rel_primes.back(), (int)vals2.cols()/rel_primes.back()/rel_primes.back())+1;
                cout << "\rComputing " << Nchi << " values of sqrt(chi)              " << flush;
                vector<CC> sqrtchi(Nchi,CC(1,0));
                for(int i=0; i<primes_list.size(); i++) {
                    int p = primes_list[i];
                    if(__gcd(p,N)>1) continue;
                    if(p>=Nchi) break;
                    CC sqrtp = sqrt(conj(fourier(chi,p%fourier.cols())));
                    for(int j=p; j<Nchi; j+=p)
                        sqrtchi[j]=sqrtchi[j/p]*sqrtp;
                }

                cout << "\rComputing matrix to diagonalize                              " << flush;
                RRMatrix Tp=RRMatrix::Zero(dim,dim);
                int l;
		int hecke_matrices_used = 0;
                for(l=0; l<primes_list.size(); l++) {
                    int p = primes_list[l];
                    if(p==2) continue;
                    if(__gcd(p,N)>1) continue;
                    if(rel_primes.back()*rel_primes.back()*p>vals2.cols()) break;
                    for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) {
                        CC x = prod3(vals2.row(0),fourier.row(chi),rel_primes[i],rel_primes[j],p,k)*RR(log(p+0.))*sqrtchi[rel_primes[i]]*sqrtchi[rel_primes[j]]*sqrtchi[p];
                        assert(x.imag()<0.00001);
                        Tp(i,j)+=x.real();
                    }
                    if(++hecke_matrices_used > 10) {
                        break;
                    }
                }

                cout << "\rConverting c to real matrix                  " << flush;
                RRMatrix cR=RRMatrix(dim,dim);
                for(int i=0; i<dim; i++) for(int j=0; j<dim; j++) {
                    CC x = c(i,j)*sqrtchi[rel_primes[i]]*sqrtchi[rel_primes[j]];
                    assert(x.imag()<0.00001);
                    cR(i,j) = x.real();
                }

                if(hecke_matrices_used>1) cout << "\rUsing " << hecke_matrices_used << " hecke matrices to be conservative         " << endl;
                cout << "\rLoading big data                  " << flush;
                CCMatrix bb(dim, M/rel_primes.back());
                for(int i=0; i<dim; i++) for(int j=0; j<M/rel_primes.back(); j++) bb(i,j)=prod2(vals2.row(0),fourier.row(chi),rel_primes[i],j,k);
                cout << "\rBasic lin alg                  " << flush;
                //CCMatrix v = (c.inverse() * Tp);

                /*bool realizable = true;
                for(int i=0; i<v.rows(); i++)
                    for(int j=0; j<v.cols(); j++)
                        if(abs((v(0,0)*v(i,j)*conj(v(0,j))*conj(v(i,0))).imag()) > RR(0.00001))
                            realizable = false;
                if(!realizable) {
                    cout << "Matrix cannot be made real" << endl;
                    cout << v << endl;
                }*/
                cout << "\rSolve for eigenvectors           " << flush;
                Eigen::GeneralizedSelfAdjointEigenSolver<RRMatrix > ces(Tp, cR);
                assert(ces.info()==Success);
                cout << "\rCompute result                        " << flush;
                RRMatrix vR = ces.eigenvectors().transpose();
                CCMatrix vC = CCMatrix(dim,dim);
                for(int i=0; i<dim; i++)
                    for(int j=0; j<dim; j++)
                        vC(i,j) = vR(i,j) * sqrtchi[rel_primes[j]];
                CCMatrix results;
                results.noalias() = vC * bb;
                adjust(results);
                //cout << results << endl;
                //cout << ces.eigenvalues() << endl;

                for(int i=0; i<results.rows(); i++) {
                    cout << "\rVerifying multiplicativity             " << flush;
                    bool good = true;
                    for(int j=0; j<results.rows(); j++) if(i!=j)
                        if(close_enough(ces.eigenvalues()(i),ces.eigenvalues()(j)))
                            good = false;
                    if(!good) {
                        cout << "\rSkipping high multiplicity eigenvalue " << ces.eigenvalues()(i) << endl;
                        number_bad++;
                        continue;
                    }
                    for(int j=2; j<results.cols(); j++) {
                      for_prime_factors(j) {
                          int j1=j, j2=1;
                          while(j1%p==0) j1/=p, j2*=p;
                          CC A = results(i,j1)*results(i,j2);
                          CC B = results(i,j);
                          if(!close_enough(A,B) && __gcd(j,N)==1) {// && __gcd(l,N)==1 ) {
                            if(true || good) {
                                cout << setprecision(0) << chi << " " << i << " " << j1 << " " << j2 << " " << A << " " << B << " " << results(i,j1) << " " << results(i,j2) << " " << results(i,j) << endl;
                                //cout << setprecision(0) << fourier.row(chi) << endl;
                                //assert(false);
                            }
                            good=false;
                          }
                          break;
                      }
                    }
                    if(!good) {
                        cout << "\rFailed multiplicativity check        " << flush;
                        number_bad++;
                    } else {
                        cout << "\rPassed multiplicativity check         " << flush;
                        number_good++;
                        char filename[100];
                        sprintf(filename, "newforms/level%d_character%d_weight%d_number%d", N, chi, k, i);
                        ofstream out(filename);
                        for(int j=1; j<results.cols(); j++)
                            out << setprecision(prec*0.28) << results(i,j).real() << " + " << results(i,j).imag() << " i" << endl;
                        out.close();
                    }
                }

                cout << "\r" << M/rel_primes.back() << " coefficients for " << number_good << " newforms (" << number_bad << " failed to be computed)" << endl;
            }

        } else {
            cout << "Calculations would overflow long, switching to gmp" << endl;
            //allTrThat12<mpz_int>(M,N,k);
            assert(false);
            abort();
        }
}
