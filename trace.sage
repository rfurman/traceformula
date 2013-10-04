# some sage code to use the trace formula
# to compute some modular forms.
#
# going to use the formula in Schoof and van der Vlugt
#
# I just want to get something working quickly, which will probably
# run slowly.

def getpoly(k):
    t,n = QQ['t,n'].gens()
    last = 1+0*t
    cur = t
    if k==0: return last
    for i in range(k-1):
        cur, last = t*cur - n*last, cur
    return cur

psi_memo = {}

def psi(N):
    if N in psi_memo:
        return psi_memo[N]
    S = N
    for (p,e) in factor(N):
        S = ZZ ( S * (p+1) / p )
    psi_memo[N] = S
    return S

# for now the character is going to be trivial...

def A1(n, N, k, chi):
    if is_square(n) and gcd(n,N)==1:
        return CC(psi(N)*(k-1)/12 * n^(k/2-1)*CC(chi(sqrt(n))))
    return 0

def A1hat(n, N, k, a):
    if is_square(n) and gcd(n,N)==1 and sqrt(n)%N==a:
        return psi(N)*(k-1)/12 * n^(k/2-1)
    return 0

hw_memo = {}

def hw(d):
    if d in hw_memo:
        return hw_memo[d] 
    D = fundamental_discriminant(d)
    f = sqrt(d/D)

    S = 1
    for p,e in factor(f):
        S = S * (1 - kronecker(D,p)/p)

    S = S * QuadraticField(D).class_number()
    if d % 4 in [2,3]:
        hw_memo[d]=0
        return 0
    if D == -3:
        S = S/3
    elif D == -4:
        S = S/2
    S = f*S
    hw_memo[d]=S
    return S

def mubad(t, f, n, N, chi):
    S = 0
    Nf = gcd(N,f)
    for x in range(N):
        if ZZ(N*Nf).divides(x*x-t*x+n) and gcd(x,N)==1:
            S+=CC(chi(x))
    S = S * psi(N)/psi(ZZ(N/gcd(N,f)))
    return S

def muhat(t, f, n, N, a):
    S = 0
    Nf = gcd(N,f)
    if (a*a-t*a+n)%(Nf*N)==0:
        return psi(N)/psi(ZZ(N/gcd(N,f)))
    else:
        return 0

def A2(n, N, poly, chi):
    S = 0
    bound = ceil(2*sqrt(n))-1
    for t in srange(-bound, bound+1):
        S2 = 0
        for f in divisors(4*n - t*t):
            if (t*t-4*n) % (f*f) == 0:
                a = ZZ((t*t - 4*n)/(f*f))
                S2 = S2 + CC(hw(a)*mubad(t, f, n, N, chi))
        S = S + S2 * poly(t,n)

    return -.5 * S

def A2hat_additive(n, N, poly, a):
    S = 0
    bound = ceil(2*sqrt(n))-1
    for t in srange(-bound, bound+1):
        S2 = 0
        for f in divisors(4*n - t*t):
            if (t*t-4*n) % (f*f) == 0:
                d = ZZ((t*t - 4*n)/(f*f))
                S2 = S2 + hw(d)*muhat(t, f, n, N, a)
        S = S + S2 * poly(t,n)

    return -S/2

cl_memo = {}

def A2hat_multiplicative(n, N, poly, y):
    S = 0
    bound = ceil(2*sqrt(n))-1
    for t in srange(-bound, bound+1):
        if not N.divides(y*y-t*y+n):
            continue
        D = fundamental_discriminant(t*t-4*n)
        if D in cl_memo:
            S2 = cl_memo[D]
        else:
            S2 = cl_memo[D] = QuadraticField(D).class_number()

        for (p,c) in factor(ZZ(sqrt((t*t-4*n)/D))):
            a = valuation(N, p)
            d = valuation(y*y-t*y+n, p)
            S3 = 0
            kDp = p-kronecker(D,p)

            if a==0:
                S3 = 1+kDp*(p^c-1)/(p-1)
            else:
                if d>=2*a and a<=c:
                    S3 = (p^a+p^(a-1)) * (1 + kDp*(p^(c-a)-1)/(p-1))
                elif c<=d-a:
                    S3 = p^c
                S3 = S3 + p^(c-1) *  kDp * max(0,min(c,a,d-a+1))
            S2 = S2 * S3

        if D == -3:
            S2 /= 3
        elif D == -4:
            S2 /= 2
        S = S + S2 * poly(t,n)

    return -S/2

def A3(n, N, k, chi):
    S1 = 0
    for d in divisors(n):
        S2 = 0
        for c in divisors(N):
            if gcd(d,c)==1 and gcd(ZZ(n/d),ZZ(N/c))==1 and gcd(ZZ(N/chi.conductor()),ZZ(n/d)-d)%gcd(c,ZZ(N/c))==0:
                y = CRT(d,ZZ(n/d),c,ZZ(N/c))
                S2 += euler_phi(gcd(c,ZZ(N/c))) * CC(chi(y))
        S1 += S2 * min(d,ZZ(n/d))^(k-1)
    return -S1/2

def A3hat(n, N, k, a):
    S1 = 0
    for d in divisors(n):
        S2 = 0
        # c divides N, divides d-a, and is a multiple of N/gcd(N,n/d-a)
        c1 = ZZ(N/(gcd(N,n/d-a)))
        c2 = gcd(N,d-a)
        #if c1.divides(c2):
        #    S2 = sigma(c2/c1,0)
        for c in divisors(gcd(N,d-a)):
            if (n/d-a)%(N/c)==0:
                g = gcd(c,N/c)
                S2 += euler_phi(g)/g
        S1 += S2 * min(d,ZZ(n/d))^(k-1)
    return -S1/2

def A4hat(n,N,k, a):
    S1 = 0
    if k==2:
        for t in divisors(n):
            if gcd(N,ZZ(n/t))==1:
                S1+=t
        S1 /= euler_phi(N)
        return S1
    else:
        return 0

def A4(n,N,k, chi):
    S1 = 0
    if k==2 and chi.is_trivial():
        for t in divisors(n):
            if gcd(N,ZZ(n/t))==1:
                S1+=t
    return S1

def TrThat(n, N, k, a):
    if n==0:
        return 0
    return A1hat(n,N,k,a)+A2hat_multiplicative(n,N,getpoly(k-2),a)+A3hat(n,N,k,a)+A4hat(n,N,k,a)


def TrT(n, N, k, chi=None):
    if chi is None:
        chi = DirichletGroup(N, base_ring=QQ)[0] # trivial character
    if n==0:
        return 0
    #if chi(-1) != (-1)^(k%2):
    #    return 0
    return A1(n,N,k,chi)+A2(n,N,getpoly(k-2),chi)+A3(n,N,k,chi)+A4(n,N,k,chi)

def allTrThat(M, N, k):
    astar = [i for i in range(N) if gcd(i,N)==1]
    return matrix( [ [ TrThat(n, N, k, a) for n in range(M) ] for a in astar ] )

def allTrT(M, N, k):
    astar = [i for i in range(N) if gcd(i,N)==1]
    #chars = [g for g in DirichletGroup(N) if g(-1)==(-1)^k]
    chars = [g for g in DirichletGroup(N)]
    fourier = matrix( [ [ char(a) for a in astar ] for char in chars ] )
    return fourier * allTrThat(M, N, k)

# coefficients of 1/zeta(s)^2
def moebius2(n):
    ret = 1
    for (p,e) in n.factor():
        if e==1:
            ret *= -2
        elif e==2:
            pass
        else:
            return 0
    return ret

# Warning! wrong answer when gcd(n,N^infty) is not square-free
def TrTnew(n, N, k, chi=None):
    if chi is None:
        Nx = 1
        chi = DirichletGroup(N, base_ring=QQ)[0]
    else:
        Nx = chi.conductor()
    if n==0:
        return 0
    sum = 0
    for d in divisors(ZZ(N/Nx)):
      if TrT(1, ZZ(N/d), k, chi) != 0:
        chi2 = chi.primitive_character().extend(N/d)
        sum += moebius2(prime_to_m_part(d,n)) * moebius(d/prime_to_m_part(d,n)) * TrT(n, ZZ(N/d), k, chi2)
        if (d*d).divides(n) and d>1:
            sum -= moebius2(prime_to_m_part(d,n)) * moebius(d/prime_to_m_part(d,n)) * TrT(n/d/d, ZZ(N/d), k, chi2) * d^3 * chi2(d)
    return sum

# Warning! wrong answer when gcd(n,N^infty) is not square-free
def allTrTnew(M, N, k):
    if chi is None:
        Nx = 1
    else:
        Nx = chi.conductor()
    if n==0:
        return 0
    sum = 0
    for d in divisors(N/Nx):
      if TrT(1, ZZ(N/d), k, chi) != 0:
        sum += moebius2(prime_to_m_part(d,n)) * moebius(d/prime_to_m_part(d,n)) * TrT(n, ZZ(N/d), k, chi)
    return sum

def test_TrT(N, k, chi=None):
    if chi is None:
        chi = DirichletGroup(N, base_ring=QQ)[0] # trivial character
    C = CuspForms(chi, k)
    if round(TrT(1,N,k,chi))==0:
        return
    print "Testing {} {} of dim {}".format(N,k,round(TrT(1,N,k,chi)))
    diffs = [(C.hecke_operator(n).trace() - TrT(n, N, k, chi))/n^((k-1)/2) for n in range(1,100)]
    if any([abs(d)>0.001 for d in diffs]):
        print "FAIL"
    else:
        print "pass"

def first_d_relatively_prime_to_n(d, n):
    gen = (i for i in xrange(1,10000000) if gcd(i,n)==1)
    return [gen.next() for i in range(d)]


def prod2(T, m, n):
    ret = 0
    for d in divisors(gcd(m,n)):
        ret += T[m*n/d/d]*d^3
    return ret

def prod3(T, m, n, p):
    ret = 0
    for d in divisors(gcd(m,n)):
        ret += prod2(T,m*n/d/d,p)*d^3
    return ret



# Given v a sum of d multiplicative sequences return the d sequences
# * v should start with a 0
# * v should be of length at least 4*p*d^2
# * so far just assumes that multiplicativity at p is enough
def multiplicative_basis(T, bad=1, p=2):
    dim = T[1].real().round()
    N = len(T)
    good = first_d_relatively_prime_to_n(dim, bad)
    c = matrix(CDF, [ [prod2(T,m,n) for n in good] for m in good] ).inverse()
    v = c * matrix(CDF, [ [prod3(T,m,n,p) for n in good] for m in good] )
    print c.inverse()
    print c
    print c.inverse()*v
    print v
    D, P = v.left_eigenmatrix()
    print D
    print P
    D, P = v.right_eigenmatrix()
    print D
    print P
    print P.transpose()
    #for n in range(1,N):
    #    T[n] *= n^1.5
    basis = P.transpose()*matrix( [ [prod2(T,m,n) for n in range(1,N/good[-1])] for m in good] )
    print P.transpose().inverse()*basis
    return [basis.row(i) / basis[i,0] for i in range(dim)]
