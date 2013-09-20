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

def psi(N):
    S = RR(1)
    for (p,e) in factor(N):
        S = S * (1.0 + 1/RR(p))
    S = S * N
    return S

# for now the character is going to be trivial...

def A1(n, N, k, chi):
    if is_square(n) and gcd(n,N)==1:
        return CC(psi(N)*(k-1)/12 * n^(k/2-1)*CC(chi(sqrt(n))))
    return 0

def hw(d):
    D = fundamental_discriminant(d)
    f = sqrt(d/D)

    S = 1
    for p,e in factor(f):
        S = S * (1 - kronecker(D,p)/p)

    S = S * QuadraticField(D).class_number()
    if d % 4 in [2,3]:
        return 0
    if D == -3:
        S = S/3
    elif D == -4:
        S = S/2
    S = f*S
    return S

def countroots(t, n, f, N):
    R1 = Integers(N*gcd(N,f))
    S = PolynomialRing(R1, 'x')
    x = S.gen()
    g = x^2 - t*n + n
    #blah ablah blahasfdadf

def mubad(t, f, n, N, chi):
    S = 0
    for x in srange(N):
        if (x*x-t*x+n)%(N*gcd(N,f))==0 and gcd(x,N)==1:
            S+=CC(chi(x))
    S = S * psi(N)/psi(N/gcd(N,f))
    return S

def mu(t, f, n, N):
    # assume N is odd (and squarefree?)
    S = 1
    a = (t^2 - 4*n)/(f*f)
    if gcd(a,N) != 1:
        return 0
    for (p,e) in factor(N * gcd(N,f)): # this probably isn't right...
        S = S * (1 + kronecker_symbol(a, p))

    S = RR(S)
    S = S * psi(N)/psi(N/gcd(N,f))
    return S

def A2(n, N, poly, chi):
    S = 0
    bound = ceil(2*sqrt(n))-1
    for t in srange(-bound, bound+1):
        S2 = 0
        for f in divisors(4*n - t*t):
            if (t*t-4*n) % (f*f) == 0:
                a = (t*t - 4*n)/(f*f)
                S2 = S2 + CC(hw(a)*mubad(t, f, n, N, chi))
        S = S + S2 * poly(t,n)

    return -.5 * S

def A3(n, N, k, chi):
    S1 = 0
    for d in divisors(n):
        S2 = 0
        for c in divisors(N):
            if gcd(d,c)==1 and gcd(n/d,N/c)==1 and gcd(N/chi.conductor(),n/d-d)%gcd(c,N/c)==0:
                y = CRT(d,ZZ(n/d),c,ZZ(N/c))
                S2 += euler_phi(gcd(c,N/c)) * CC(chi(y))
        S1 += S2 * min(d,n/d)^(k-1)
    return -S1/2

def A4(n,N,k, chi):
    S1 = 0
    if k==2 and chi.is_trivial():
        for t in divisors(n):
            if gcd(N,n/t)==1:
                S1+=t
    return S1

def TrT(n, N, k, chi=None):
    if chi is None:
        chi = DirichletGroup(N, base_ring=QQ)[0] # trivial character
    if n==0:
        return 0
    if chi(-1) != (-1)^(k%2):
        return 0
    return A1(n,N,k,chi)+A2(n,N,getpoly(k-2),chi)+A3(n,N,k,chi)+A4(n,N,k,chi)

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
    else:
        Nx = chi.conductor()
    if n==0:
        return 0
    sum = 0
    for d in divisors(N/Nx):
        sum += moebius2(prime_to_m_part(d,n)) * moebius(d/prime_to_m_part(d,n)) * TrT(n, N/d, k, chi)
    return sum

def first_d_relatively_prime_to_n(d, n):
    gen = (i for i in xrange(10000000) if gcd(i,n)==1)
    return [gen.next() for i in range(d)]

# Given v a sum of d multiplicative sequences return the d sequences
# * v should start with a 0
# * v should be of length at least 4*p*d^2
# * so far just assumes that multiplicativity at p is enough
def multiplicative_basis(T, dim, bad=1, p=2):
    good = first_d_relatively_prime_to_n(dim, bad*p)
    N = len(T)
    c = matrix(CDF, [ [T[m*n] for n in good] for m in good] ).inverse()
    v = c * matrix(CDF, [ [T[m*n*p] for n in good] for m in good] )
    D, P = v.right_eigenmatrix()
    basis = P.transpose() * matrix( [ [T[m*n] for n in range(1,N/good[-1])] for m in good] )
    basis = diagonal_matrix([1/basis[i,i] for i in range(dim)]) * basis
    return basis.rows()
