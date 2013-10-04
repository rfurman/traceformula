#define for_prime_factors(N) for(int p=2; (p*p<=(N) || (p=(N))) && p>1; p++) if((N)%p==0)
int rel_prime[1000];
vector<int> primes_list;
int mobius[500000];
int mobius2[500000];

long pow(long x, int y) {
    if(y==0) return 1;
    long ret = pow(x,y/2);
    ret *= ret;
    return (y%2?ret*x:ret);
}


int pow(int x, int y) {
    if(y==0) return 1;
    long ret = pow(x,y/2);
    ret *= ret;
    return (y%2?ret*x:ret);
}

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


int classnumbers[10000000];
int classnumber(int D) {
    int& ret = classnumbers[-D];
    if(ret) return ret;
    pari_sp av = avma; ret=itos(classno(stoi(D))); avma = av;
    return ret;
}
