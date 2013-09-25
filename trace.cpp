#include <pari/pari.h>
#include <algorithm>
#include <iostream>
#include <cassert>
typedef long long i64;
using namespace std;

#define POLY(I, EXPR) int poly##I(int t, int n) { return EXPR; }
#define for_prime_factors(N) for(int p=2; (p*p<=(N) || (p=(N))) && p>1; p++) if((N)%p==0)

POLY(0, 1);
POLY(1, t);
POLY(2, t*t-n);
POLY(3, t*t*t-2*t*n);
POLY(4, t*t*t*t - 3*t*t*n + n*n);
POLY(5, t*t*t*t*t - 4*t*t*t*n + 3*t*n*n);
int (*poly[6]) (int t, int n) = { poly0, poly1, poly2, poly3, poly4, poly5 };

int psi(int N) {
    int ret=N;
    for_prime_factors(N) {
        ret = ret * (p+1)/p;
        do { N/=p; } while(N%p==0);
    }
    return ret;
}

int phi(int N) {
    int ret=N;
    for_prime_factors(N) {
        ret = ret * (p-1)/p;
        do { N/=p; } while(N%p==0);
    }
    return ret;
}

int number_of_divisors(int N) {
    int ret=1;
    for_prime_factors(N) {
        int e=0;
        do { N/=p; e++; } while(N%p==0);
        ret *= e+1;
    }
    return ret;
}

int sqfree(int N) {
    int ret=1;
    for_prime_factors(N) {
        int e=0;
        do { N/=p; e++; } while(N%p==0);
        if(e%2==1) ret *= p;
    }
    return ret;
}

int valuation(int N, int p) {
    int ret=0;
    while(N%p==0) N/=p, ret++;
    return ret;
}

int evalpoly(int k, int t, int n) {
    int val[2];
    val[0]=1;
    val[1]=t;
    if(k<2) return val[k];
    for(int i=2; i<=k; i++) {
        val[i%2] = t*val[(i+1)%2] - n*val[i%2];
    }
    return val[k%2];
}

int classnumbers[100000];
int classnumber(int D) {
    int& ret = classnumbers[-D];
    if(ret) return ret;
    pari_sp av = avma; ret=itos(classno(stoi(D))); avma = av;
    return ret;
}

void allTrThat12(int M, int N, int k) {
    i64 vals[N][M];
    memset(vals,0,sizeof(vals));
    int phiN = phi(N);

    // A1
    int psiN = psi(N);
    if(true) for(int n=1; n*n<M; n++) if(__gcd(n,N)==1) {
        vals[n%N][n*n] += phiN * psiN * (k-1) * pow(n, k-2);
    }

    // A2
    int bound=sqrt(4*M)+3;
    if(true) for(int t=-bound; t<=bound; t++) {
        for(int y=0; y<N; y++) {
            int n = ((t*y-y*y)%N+N)%N;
            n += (t*t/4)/N*N-N;
            while(t*t>=4*n) n+=N;
            for(; n<M; n+=N) {
                i64 D = -sqfree(4*n-t*t);
                if(D%4==-2 || D%4==-1) D*=4;
                int S2 = classnumber(D);
                int X = round(sqrt((t*t-4*n)/D));
                for_prime_factors(X) {
                    int c=0;
                    do { X/=p; c++; } while(X%p==0);
                    int a = valuation(N, p);
                    int d = valuation(y*y-t*y+n, p);
                    int S3 = 0;
                    int kDp = p-kross(D,p);

                    if(a==0) {
                        S3 = 1+kDp*(pow(p,c)-1)/(p-1);
                    } else {
                        if(d>=2*a && a<=c)
                            S3 = pow(p,a-1)*(p+1)*(1+kDp*(pow(p,c-a)-1)/(p-1));
                        else if(c<=d-a)
                            S3 = pow(p,c);
                        S3 += pow(p,c-1) * kDp * max(0,min(min(c,a),d-a+1));
                    }
                    S2 *= S3;
                }
                S2 *= 6 * phiN;
                if(D==-3) S2/=3;
                else if(D==-4) S2/=2;
                vals[y][n] -= S2*evalpoly(k-2, t, n);
                //if(t>0) vals[N-1-y][n] -= S2*evalpoly(k-2, -t, n);
            }
        }
    }

    // A3
    if(true) for(int d=1; d<M; d++) {
        for(int n=d; n<M; n+=d) {
            for(int a=0; a<N; a++) {
                i64 c1=N/__gcd(N,abs(n/d-a));
                i64 c2 = __gcd(N,abs(d-a));
                if(c2%c1==0) {
                    vals[a][n] -= phiN*6*number_of_divisors(c2/c1)*pow(min(d,n/d),k-1);
                }
            }
        }
    }

    // A4
    if(true) if(k==2) {
        for(int t=1; t<M; t++) {
            for(int a=0; a<N; a++) {
                for(int n=t; n<M; n+=t) {
                    vals[a][n] += 12*t;
                }
            }
        }
    }

    for(int a=0; a<N; a++) {
        for(int n=1; n<M; n++) {
            cout << vals[a][n] << " ";
        }
        cout << endl;
    }
}

int main(void) {
        pari_init(8000000, 500000);

        int M=100;

        allTrThat12(40, 9, 4);
}
