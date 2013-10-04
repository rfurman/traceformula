// Evaluate the polynomials from the trace formula
// If t^2,n <= M then bounded by Fibonacci_k * M^(k/2)a < (3M)^(k/2)
template<class i64> i64 evalpoly(int k, long t, long n) {
    switch(k) {
    case 0: return (i64)1;
    case 1: return (i64)t;
    case 2: return (i64)t*t-n;
    case 3: return (i64)t*t*t-2*t*n;
    case 4: return (i64)t*t*t*t - 3*t*t*n + n*n;
    case 5: return (i64)t*t*t*t*t - 4*t*t*t*n + 3*t*n*n;
    case 10:
        i64 tt = t*t;
        return - n*n*n*n*n + tt * (15*n*n*n*n + tt * (-35*n*n*n + tt * (28*n*n + tt * (-9*n + tt))));
    }

    i64 val[2];
    val[0]=t*t*t*t - 3*t*t*n + n*n;
    val[1]=t*t*t*t*t - 4*t*t*t*n + 3*t*n*n;
    for(int i=6; i<=k; i++) {
        val[i%2] = t*val[(i+1)%2] - n*val[i%2];
    }
    return val[k%2];
}

double theoretical_TrT_bound(int M, int N, int k) {
    return pow(double(12*M),(k-1)/2.)*N*N*log(1+M);
}

template<class T> void serialize_matrix_to_file(char *filename, const Matrix<T,Dynamic,Dynamic,ColMajor>& mat) {
    ofstream out(filename, ios::out|ios::binary);
    int bytes_to_write = mat.rows()*mat.cols()*sizeof(T);
    cout << "\rSaving to file " << filename << "                     " << flush;
    out.write((char *)mat.data(), bytes_to_write);
    out.close();
}

template<class T> bool try_deserialize_matrix_from_file(char *filename, Matrix<T,Dynamic,Dynamic,ColMajor>& mat) {
    ifstream in(filename, ios::in|ios::binary|ios::ate);
    if(in.is_open()) {
        int size = in.tellg();
        int bytes_to_read = mat.rows()*mat.cols()*sizeof(T);
        if(size >= bytes_to_read) {
            cout << "\rLoading from file " << filename << "                   " << endl;
            in.seekg(0, ios::beg);
            in.read((char *)mat.data(), bytes_to_read);
            in.close();
            return true;
        } else {
            cout << "\rFile " << filename << " too small                " << flush;
            return false;
        }
    }
    return false;
}

// Compute 12 times the fourier transform of the traces
template<class i64> ZZMatrix allTrThat12(const int M, const int N, int k) {

    const int phiN = phi(N);
    //vector<vector<i64> > vals(phiN,vector<i64>(M));
    ZZMatrix vals = ZZMatrix::Zero(N,M);

    char filename[100];
    sprintf(filename, "trace_hats/level%d_weight%d", N, k);
    if(try_deserialize_matrix_from_file(filename,vals)) return vals;

    for(int i=0, j=0; i<N; i++) if(__gcd(i,N)==1) rel_prime[j++]=i;
    rel_prime[phiN]=0;

    int msk = 15;

    // A1
    // Bounded by N^2 k M^(k/2-1)
    int psiN = psi(N);
    if(msk&1) for(int i=0; i<phiN; i++) for(int n=rel_prime[i]; n*n<M; n+=N) {
        vals.coeffRef(rel_prime[i],n*n) += (phiN * psiN * (k-1)) * pow((i64)n, k-2);
    }

    // A2
    // Bounded by (12M)^(k/2-1) * 6N^2 * h(4M) * sqrt(M)/N < (12M)^(k/2) * N * log log M
    int bound=sqrt(4*M)+1;
    if(msk&2) for(int t=0; t<=bound; t++) {
        for(int i=0; i<phiN; i++) {
            int y=rel_prime[i];
            int n = ((t*y-y*y)%N+N)%N;
            n += (t*t/4)/N*N-N;
            while(t*t>=4*n) n+=N;
            for(; n<M; n+=N) {
                int D = -sqfree(4*n-t*t);
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
                        if(d>=2*a && a<=c) S3 = pow(p,a-1)*(p+1)*(1+kDp*(pow(p,c-a)-1)/(p-1));
                        else if(c<=d-a) S3 = pow(p,c);
                        S3 += pow(p,c-1) * kDp * max(0,min(min(c,a),d-a+1));
                    }
                    S2 *= S3;
                }
                S2 *= 6 * phiN;
                if(D==-3) S2/=3;
                else if(D==-4) S2/=2;
                vals.coeffRef(y,n) -= S2*evalpoly<i64>(k-2, t, n);
                if(t>0) vals.coeffRef((N-y)%N,n) -= S2*evalpoly<i64>(k-2, -t, n);
            }
        }
    }

    // A3
    // Bounded by 12*M^(k/2+1/2)*log(M)*log(N)
    if(msk&4) for(int d=1; d<M; d++) {
        for(int n=d; n<min(M,d*d+1); n+=d) {
            i64 p = pow((i64)min(d,n/d),k-1);
            int S2=0;
            for(int i=0; i<phiN; i++) {
                int a=rel_prime[i];
                int c1=N/__gcd(N,abs(n/d-a));
                int c2 = __gcd(N,abs(d-a));
                if(c2%c1==0) {
                    for(int c=c1; c<=c2; c+=c1) if(c2%c==0) {
                        int g=__gcd(c,N/c);
                        vals.coeffRef(a,n) -= (d*d==n?6:12)*p*phi(g)*phi(N/g);
                    }
                }
            }
        }
    }

    // A4
    // Bounded by N^2*log(N)
    if(msk&8) if(k==2) {
        for(int t=1; t<M; t++) {
            for(int i=0; i<phiN; i++) {
                int a = rel_prime[i];
                for(int n=t; n<M; n+=t) {
                    vals.coeffRef(a,n) += 12*t;
                }
            }
        }
    }

    double theoretical_upper_bound = theoretical_TrT_bound(M,N,k);
    double actual_upper_bound = 0;

    for(int a=0; a<phiN; a++) {
        for(int n=1; n<M; n++) {
            //double x = static_cast<double>(vals.coeffRef(a,n));
            double x = to_double(vals.coeffRef(rel_prime[a],n));
            if(x>actual_upper_bound) {
                actual_upper_bound = x;
            }
        }
    }
    cout << "\r" << theoretical_upper_bound << " theoretical bound vs " << actual_upper_bound << " actual        " << endl;
    assert(theoretical_upper_bound > actual_upper_bound);

    serialize_matrix_to_file(filename,vals);

    return vals;
}

// Compute 12 times the fourier transform of the traces, sieved for newforms
RRMatrix allTrThat12new(int M, int N, int k) {
    RRMatrix vals = RRMatrix::Zero(N,M);
    for(int d=1; d<=N; d++) if(N%d==0) {
        cout << "\rComputing forms of level " << N/d << "             " << endl;
        ZZMatrix vals2 = allTrThat12<long>(M,N/d,k);
        for(int y=0; y<N; y++) if(true||__gcd(N,y)==1) for(int n=0; n<M; n++) {
            int prime_to_n_part = d;
            while(__gcd(prime_to_n_part,n)>1) prime_to_n_part /= __gcd(prime_to_n_part,n);
            vals(y,n) += mobius2[prime_to_n_part] * mobius[d/prime_to_n_part] * vals2(y%(N/d),n);
            if(d>1 && n%(d*d)==0 && N%(d*d)!=0) {
                int dinv=0;
                for(int i=0; i<N/d; i++) if((d*i)%(N/d)==1) dinv=i;
                //if(N/d>1) assert(dinv);
                int dk=1; // d^(k-1)
                for(int i=0; i<k-1; i++) dk*=d;
                vals(y,n) -= mobius2[prime_to_n_part] * mobius[d/prime_to_n_part] * vals2((y*dinv)%(N/d),n/d/d) * dk;
            }
        }
    }

    return vals;
}

