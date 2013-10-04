#include "slint.h"

using namespace std;

inline CC e(RR z) {
    CC twopii(0,2*mpfr::const_pi());
    return exp(twopii * z);
}

class DirichletGroup;


class DirichletCharacter {
public:
    int m;         // the label
    DirichletGroup * parent;

    DirichletCharacter() {}
    DirichletCharacter(DirichletGroup * parent_, int m_);

    ~DirichletCharacter() {}
    int exponent(int n);
    CC max(int* index);
    CC value(int n);
    bool is_primitive();
    bool is_primitive_at_two();
    bool is_primitive_at_odd_part();
    bool is_even();

    CC gauss_sum();
};


class DirichletGroup {
    //
    //
public:
    int q;         // the modulus
    int q_odd;     // the odd part of the modulus
    int q_even;

    int k;          // the number of odd prime factors of q
 
    std::vector<int> * primes;
    std::vector<int> * exponents;
    std::vector<int> * generators;

    int ** A;      // exponent vectors:
                    //      for each m coprime to q we store an array
                    //      with the property that
                    //
                    //          m == g[j]**A[m][k] mod p[j]**e[j]
                    //
                    //      where g[j] is a primitive root mod p[j]**e[j],
                    //      and p[j] is the j-th prime factor of q.
                    //      (We don't actually store g[k], p[j], or e[j] anywhere.)
    
    int * B;

    int * PHI;     // PHI[j] = phi(q/phi(p[j]**e[j])). This will make it easier
                    // to compute the characters.

    int phi_q_odd;
    int phi_q;     // phi(q)

    CC * zeta_powers_odd;     // zeta_powers[n] == e(n/phi(q))
    CC * zeta_powers_even;     // zeta_powers[n] == e(n/phi(q))
    
    DirichletGroup(int q_) : q(q_) {
        q_even = 1;
        q_odd = q;
        while(q_odd % 2 == 0) {
            q_odd /= 2;
            q_even *= 2;
        }

        phi_q = euler_phi(q);

        if(q_odd > 1) {
            primes = new std::vector<int>();
            exponents = new std::vector<int>();
            generators = new std::vector<int>();

            factors(q_odd, primes, exponents);
            k = primes->size();
            phi_q_odd = euler_phi(q_odd);
            
            PHI = new int[k];
            A = new int*[q_odd];
            for(int n = 0; n < q_odd; n++) {
                A[n] = new int[k];
            }
            zeta_powers_odd = new CC[phi_q_odd];

            for(int j = 0; j < k; j++) {
                int x = pow(primes->at(j), exponents->at(j));
                int g = primitive_root(x);
                generators->push_back(g);
                int phi = (pow(primes->at(j), exponents->at(j) - 1) * (primes->at(j) - 1));
                PHI[j] = phi_q_odd/phi;
                int a = 1;
                for(int l = 0; l < phi; l++) {
                    for(int m = a; m < q_odd; m += x) {
                        A[m][j] = l;
                    }
                    a = (a * g) % x;
                }
            }
        
            //
            // for each m, 0 <= m < q, if (m,q) > 1, store
            // this as a flag in A[m][0]
            //
            for(int m = 0; m < q_odd; m++) {
                if(GCD(m,q_odd) > 1) {
                    A[m][0] = -1;
                }
            }

            for(int n = 0; n < phi_q_odd; n++) {
                zeta_powers_odd[n] = e(n/(RR)phi_q_odd);
            }
        } // end of initialization of everything having to do
          // with q_odd
        
        if(q_even > 4) {
            B = new int[q_even];
            zeta_powers_even = new CC[q_even/4];
            for(int n = 0; n < q_even/4; n++) {
                zeta_powers_even[n] = e(4*n/RR(q_even));
            }
            int pow_three = 1;
            for(int e = 0; e < q_even/4; e++) {
                B[pow_three] = e;
                B[pow_three - 1] = 1;
                B[q_even - pow_three] = e;
                B[q_even - pow_three - 1] = -1;
                pow_three *= 3;
                pow_three %= q_even;
            }
        }

    }

    ~DirichletGroup() {
        if(q_odd > 1) {
            delete [] zeta_powers_odd;
            for(int n = 0; n < q_odd; n++) {
                delete [] A[n];
            }
            delete [] A;
            delete [] PHI;
            delete primes;
            delete exponents;
            delete generators;
        }
        if(q_even > 4) {
            delete [] B;
            delete [] zeta_powers_even;
        }
    }

    int chi_odd_exponent(int m, int n) {
        int x = 0;
        for(int j = 0; j < k; j++) {
            x += A[m][j]*A[n][j]*PHI[j];
            //x = x % phi_q_odd;
            if(x > 1000000000)
                x = x % phi_q_odd;
        }
        if(x >= phi_q_odd) 
            x = x % phi_q_odd;

        return x;
    }
    
    int chi_even_exponent(int m, int n) {
        int x = B[m]*B[n];
        if(B[m-1] == -1 && B[n-1] == -1)
            x += q_even/8;
        return x % (q_even/4);
    }

    CC chi(int m, int n) {
        CC even_part = (RR)1;
        CC odd_part = (RR)1;
        if(q_even > 1) {
            if(m % 2 == 0 || n % 2 == 0) {
                return (RR)0;
            }
            else if(q_even == 2) even_part = 1;
            else if(q_even == 4) {
                if(m % 4 == 3 && n % 4 == 3) even_part = -1;
                else even_part = 1;
            }
            else {
                even_part = zeta_powers_even[chi_even_exponent(m % q_even, n % q_even)];
            }
        }
        if(m >= q_odd)
            m %= q_odd;
        if(n >= q_odd);
            n %= q_odd;
        if(q_odd > 1) {
            if(A[m][0] == -1 || A[n][0] == -1)
                return (RR)0;
            else
                odd_part = zeta_powers_odd[chi_odd_exponent(m, n)];
            if(q_even == 1)
                return odd_part;
        }
        return odd_part * even_part;
    }

    //DirichletCharacter character(int m) {
    //    return DirichletCharacter(q, m, k, A, PHI, phi_q, zeta_powers);
    //}
    DirichletCharacter character(int m) {
        return DirichletCharacter(this, m);
    }

    //DirichletCharacter * character(int n) {
        //
        // return the character chi for which
        // chi(a) == e(n/(q-1))
        //

    //    return new DirichletCharacter(q, n, zeta_powers, power_table);
    //}


    int exponent(int n, int m) {
        //
        //Return the number a such that chi(n, m) = e(a/phi(q)).
        //
        int e;
        int odd_exponent;
        int even_exponent;

        if(q_odd > 1) {
            odd_exponent = chi_odd_exponent(n % q_odd, m % q_odd);
        }
        else {
            odd_exponent = 0;
        }


        if(q_even > 4) {
            even_exponent = chi_even_exponent(n % q_even, m % q_even);
            even_exponent *= 2; // the function just above computes the exponent of
                                // e(1/ (q_even/4) ), but we want the exponent of
                                // e(1/phi(q_even)) = e(1/(q_even/2))
        }
        else if(q_even == 4) {
            if ((n % q_even) == 3 && (m % q_even) == 3) {
                even_exponent = 1;
            }
            else {            
                even_exponent = 0;
            }
        }
        else {
            even_exponent = 0;
        }
        
        if(q_even == 1) { // special case because phi(1) != 1/2.
            e = odd_exponent;
        }
        else {
            e = odd_exponent * q_even/2 + even_exponent * phi_q_odd;
        }
    
        // we now have the value of chi(m) as e(exponent/phi(q))

        // it could be equal to phi(q), though, and in that case we
        // want it to be zero...
        if(e == phi_q){
            e -= phi_q;
        }

        return e;
    }
};




DirichletCharacter::DirichletCharacter(DirichletGroup * parent_, int m_) : parent(parent_), m(m_) {
        //
        // This doesn't actually do anything. Computing characters
        // using DirichletCharacter is not going to be any faster
        // than using DirichletGroup::chi(), at least for now.
        //
}

int DirichletCharacter::exponent(int n) {
    int x = 0;
    for(int j = 0; j < parent->k; j++) {
        x += parent->A[m][j]*parent->A[n][j]*parent->PHI[j];
        if(x > 1000000000)
            x = x % parent->phi_q;
    }
    if(x >= parent->phi_q) {
        x = x % parent->phi_q;
    }
    return x;
}

CC DirichletCharacter::value(int n) {
    return parent->chi(m,n);
    //if(parent->A[m][0] == -1 || parent->A[n][0] == -1)
    //    return 0;
    //return parent->zeta_powers[exponent(n)];
}

bool DirichletCharacter::is_primitive() {
    // return whether or not this character is primitive

    return is_primitive_at_two() && is_primitive_at_odd_part();
}

bool DirichletCharacter::is_primitive_at_odd_part() {
    // this is computed one prime at a time, and the
    // character will be primitive if and only if it
    // is primitive at every prime.
    if(parent->q_odd == 1) return true;
    else {
        int n = m % parent->q_odd;
        if(parent->A[n][0] == -1) return false;
        for(int j = 0; j < parent->k; j++) {
            if(parent->A[n][j] % parent->primes->at(j) == 0)
                return false;
        }
        return true;
    }
}

bool DirichletCharacter::is_primitive_at_two() {
    int q_even = parent->q_even;
    int * B = parent->B;
    int n = m % q_even;
    if(q_even == 1) return true;
    else if(q_even == 2) return false;
    else if(q_even == 4) return n == 3;
    else if(q_even == 8) {
        if(B[n] % 2 == 1 && B[n-1] == 1) return true;
        else if(B[n] % 2 == 0 && B[n-1] == -1) return true;
        else return false;
    }
    else return B[n] % 2 == 1;
}

bool DirichletCharacter::is_even() {
    // return whether or not this character is even
    //

    //
    // We just figure out if the character is even by evaluating it at q-1.
    // The evaluation isn't going to be exact, but since the number is just
    // +-1 we can just check which of these it is closest to.
    //
    
    return abs(value(parent->q - 1) - (RR)1) < .5;
}

// Return a matrix of the character values for the characters with the same parity as k
CCMatrix matrix_of_characters_by_parity(int q, int k) {
    int nchi = 0;
    DirichletGroup dg(q);
    CCMatrix mat(q,q);
    for(int i=1; i<=dg.phi_q; i++) {
        if(dg.character(i).is_even() ^ (k%2)) {
            for(int j=0; j<q; j++) {
                mat(nchi,j) = dg.chi(i,j);
            }
            nchi++;
        }
    }
    mat.conservativeResize(nchi,q);
    return mat;
}
