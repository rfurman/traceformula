# some sage code to use the trace formula
# to compute some modular forms.
#
# going to use the formula in Schoof and van der Vlugt
#
# I just want to get something working quickly, which will probably
# run slowly.

R = RealField(53)

def getpoly(k):
    t,n = QQ['t,n'].gens()
    last = 1+0*t
    cur = t
    if k==0: return last
    for i in range(k-1):
        cur, last = t*cur - n*last, cur
    return cur


def psi(N):
    S = R(1)
    for (p,e) in factor(N):
        S = S * (1.0 + 1/R(p))
    S = S * N
    return S

# for now the weight is going to be 2 and the character
# is going to be trivial...

def A1(n, N, k):
    if is_square(n) and gcd(n,N)==1:
        return psi(N)*(k-1)/12 * n^(k/2-1)
    return 0

def hw0(d):
    if d == -3:
        return R(1/3)
    elif d == -4:
        return R(1/2)
    elif d % 4 in [2,3]:
        return 0
    else:
        #return QuadraticField(d).class_number()
        #dd = squarefree_part(d)
        #return sage.quadratic_forms.special_values.quadratic_L_function__exact(1,d)*sqrt(abs(d))/pi
        return len(sage.quadratic_forms.binary_qf.BinaryQF_reduced_representatives(d,True))

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


def H(d):
    if d == 0:
        return R(-1/12)
    if d % 4 in [1,2]:
        return 0
    if d < 0:
        return 0
    ret = 0
    for f in divisors(d):
        if d%(f*f)==0:
            ret += hw(-d/(f*f))
    return ret
    

def H1(d):
    if d == 0:
        return R(-1/12)
    if d % 4 in [1,2]:
        return 0
    if d < 0:
        return 0
    ret = 0
    for f in divisors(d):
        if d%(f*f)==0:
            ret += f*hw(-d/(f*f))
    return ret
    

def countroots(t, n, f, N):
    R1 = Integers(N*gcd(N,f))
    S = PolynomialRing(R1, 'x')
    x = S.gen()
    g = x^2 - t*n + n
    #blah ablah blahasfdadf

def mubad(t, f, n, N):
    S = 0
    for x in range(N):
        if (x*x-t*x+n)%(N*gcd(N,f))==0 and gcd(x,N)==1:
            S+=1
    S = S * psi(N)/psi(N/gcd(N,f))
    return S

def mu(t, f, n, N):
    # assume N is odd (and squarefree?)
    S = 1
    a = (t^2 - 4*n)/(f*f)
    if gcd(a,N) is not 1:
        return 0
    for (p,e) in factor(N * gcd(N,f)): # this probably isn't right...
        S = S * (1 + kronecker_symbol(a, p))

    S = R(S)
    S = S * psi(N)/psi(N/gcd(N,f))
    return S

def A2(n, N, poly):
    S = 0
    bound = ceil(2*sqrt(n))-1
    for t in range(-bound, bound+1):
        S2 = 0
        for f in divisors(4*n - t*t):
            if (t*t-4*n) % (f*f) == 0:
                a = (t*t - 4*n)/(f*f)
                S2 = S2 + RR(hw(a)*mubad(t, f, n, N))
        S = S + S2 * poly(t,n)

    return -.5 * S

def A3(n, N, k):
    S1 = 0
    for d in divisors(n):
        S2 = 0
        for c in divisors(N):
             if gcd(d,c)==1 and gcd(n/d,N/c)==1 and gcd(N,n/d-d)%gcd(c,N/c)==0:
                 S2 += euler_phi(gcd(c,N/c))
        #    if gcd(N/1,n/d-d) % gcd(c,N/c) == 0:
        #        S2 += euler_phi(gcd(c,N/c))
        S1 += S2 * min(d,n/d)^(k-1)
    return -S1/2

def A4(n,N,k):
    S1 = 0
    if k==2: # and chi==1
        for t in divisors(n):
            if gcd(N,n/t)==1:
                S1+=t
    return S1

def TrT(n, N, k):
    return A1(n,N,k)+A2(n,N,getpoly(k-2))+A3(n,N,k)+A4(n,N,k)



def TrT1(n, N):
    S2 = 0
    bound = ceil(2*sqrt(n))
    for m in range(-bound, bound+1):
        S2 += H(4*n-m*m) * (m*m-n)
    S2 *= -1/2

    S4 = 0
    for d in divisors(n):
        S4 += min(d,n/d)^(k-1)
    S4 *= -1/2

    return S2+S4
