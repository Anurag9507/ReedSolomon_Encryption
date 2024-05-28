# IMT2022103 Anurag Ramaswamy
# IMT2022090 R Lakshman

import gmpy2
from gmpy2 import mpz
import random

# All the helper functions used in our implementation

def rand(a, b):
    return mpz(a + random.random() * (b - a)) 
# Returns a value from a to b-1 (including both)

def fast_pow(a,b): # Function that computes a^b using binary exponentiation
    ans = 1
    while(b>0):
        if(mod(b,2)!=0):
            ans = ans * a
        a = a * a
        b = b >> 1
    return ans

def mod(a, b):
    ans = a - b * (a // b) # Remainder when a is divided by b
    if(ans < 0): # Adding or subtracting b to ensure ans is positive 
        if(b > 0):
            ans = ans + b
        else:
            ans = ans - b 
    return ans

def fast_pow_mod(a,b,n): # Function that returns a^b mod n (combining above 2 functions) 
    a = mod(a,n)
    ans=1
    while(b>0):
        if(mod(b,2)!=0):
            ans = mod(ans * a, n)
        a = mod(a * a, n)
        b=b>>1
    return ans

def mod_inverse(a, m):
    gcd, x, y = binaryEGCD(a, m)
    if gcd != 1: # Then mod inverse does not exist
        raise ValueError("The modular inverse does not exist.")
    else:
        return mod(x,m)  # This will be the modular inverse
    
def get_k_primes(k): # Returns the first k odd primes
    primes = []
    n = 7993632487
    while len(primes) < k:
        if miller_rabin(n,10):
            primes.append(n)
        n += 2
    return primes

def binaryEGCD(a,b): # Binary EGCD algorithm
    r=a
    r1=b
    e=0
    while(mod(r,2) == 0 and mod(r1,2) == 0):
        r=r >> 1
        r1=r1 >> 1
        e=e+1
    a1,b1,s,t,s1,t1=r,r1,1,0,0,1
    while True:
        while(mod(r,2) == 0):
            r = r >> 1
            if(mod(s,2) == 0 and mod(t,2) == 0):
                s,t=s >> 1, t >> 1
            else:
                s,t=(s+b1) >> 1,(t-a1) >> 1
        while(mod(r1,2) == 0):
            r1 = r1 >> 1
            if(mod(s1,2) == 0 and mod(t1,2) == 0):
                s1,t1=s1 >> 1,t1 >> 1
            else:
                s1,t1=(s1+b1) >> 1,(t1-a1) >> 1
        if(r1 < r):
            r, s, t, r1, s1, t1 = r1, s1, t1, r, s, t
        r1 = r1 - r
        s1 = s1 - s
        t1 = t1 - t
        if(r1 == 0):
            break
    gcd = fast_pow(2,e) * r
    return gcd, s, t

def crt(n_list, a_list): # Chinese remainder theorem (taking a list of n's and a's)
    x = 0
    N = mpz(1)
    for n_i in n_list:
        N *= n_i
    length=len(n_list)
    nstar=[gmpy2.mpz(1) for i in range(length)]
    e=[gmpy2.mpz(1) for i in range(length)]
    for i in range(length):
        nstar[i]= N // n_list[i]
        b = mod(nstar[i],n_list[i])
        t = mod_inverse(b,n_list[i])
        e[i] = nstar[i] * t
    ans = gmpy2.mpz(0)
    for i in range(length):
        ans=mod(ans + a_list[i]*e[i],N)
    return ans
        

def miller_rabin(n,rounds):
#MR2 algorithm, primality testing algorithm with a threshold value for tests while running
    if(n == 1):
        return False
    if(n == 2 or n == 3):
        return True
    if(mod(n,2) == 0 ):
        return False
    t,h=n-1,1
    while(mod(t,2) == 0):
        t = t >> 1
        h = h + 1
    while(rounds!=0):
        alpha=rand(2,n-1)
        beta=fast_pow_mod(alpha,t,n)
        for j in range(h):
            y = fast_pow_mod(beta,2,n)
            if(y == 1 and beta != 1 and beta != n-1):
                return False
            beta = y
        if(y != 1):
            return False
        rounds=rounds-1
    return True
    
def c(x): 
    return mpz(x) # Convert a number into mpz object


def extended_euclidean(a, b): 
# Extended euclidean algorithm while storing all ri s, si s, ti s and returning them 
    r, s, t = a, c(1), c(0)
    r0, s0, t0 = b, c(0), c(1)
    # lists to store values of r, s, t at each iteration
    r_values, s_values, t_values = [r], [s], [t]
    while r0 != 0:
        q = r // r0
        r, r0 = r0, mod(r,r0)
        s, s0 = s0, s - s0 * q
        t, t0 = t0, t - t0 * q   
        r_values.append(r)
        s_values.append(s)
        t_values.append(t)
    return r, s, t, r_values, s_values, t_values # r is gcd   