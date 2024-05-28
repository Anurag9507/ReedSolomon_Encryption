# IMT2022103 Anurag Ramaswamy
# IMT2022090 R Lakshman

from helper import *
import math
import gmpy2
from gmpy2 import mpz

M = mpz("1"+"0"*1000) # Upper bound on the message we are sending
mu = 0.3 # Fraction of values that can get corrupted
k = 1000 # Number of residues transmitted 

primes=[] # To store k primes 

def GlobalSetup(mu, M):
    global primes
    global k
    primes=get_k_primes(k)

def ReedSolomonSend(a):
    mods=[]
    for i in primes:
        mods.append(mod(a,i))  #mods list is a mod ni (here ni list is our primes list)
    b=Transmit(mods)
    return b

def Transmit(mods):
    global l
    global k
    l = rand(1,mu*k+1)    # L is how many values are going to get corrupted 
    l_copy=l
    corrupted_indexes=[] # Randomly choosing l corrupted indices
    while(l_copy > 0):
        random_index=rand(0,k)
        if(random_index not in corrupted_indexes):
            corrupted_indexes.append(random_index)
            l_copy = l_copy  - 1
    b=mods.copy()
    for i in range (0,k):
        if(mpz(i) in corrupted_indexes):
            p_i=primes[i]
            rand_int=rand(0,p_i)
            while(rand_int == mods[i]):
                rand_int=rand(0,p_i)
            b[i]=rand_int
    return b # This will be the corrupted list of residues recieved at the other end 

def ReedSolomonReceive(b): # implemeting Reed Solomon algorithm
    global primes
    B=crt(primes,b) # N is product of all selected primes 
    
    n=1
    for x in primes:
        n*=x
    r, s, t, r_values, s_values, t_values = extended_euclidean(n,B) 
    l_largest_primes=[] 
    primes.sort(reverse=True)
    l_largest_primes=primes[0:l]
    P=1
    for pi in l_largest_primes: 
        P=P*pi  # Product of l largest primes
    global M
    rstar = M*P                
    tstar = P
    ind = 0
    while(r_values[ind]>rstar):      # First index j such that r[j] > rstar 
        ind = ind+1   
                          
    r1 = r_values[ind] 
    s1 = s_values[ind] 
    t1 = t_values[ind]  
    if mod(r1,t1)==0: 
        return r1 // t1; 
    else: 
        return 'Error'
      
GlobalSetup(mu,M)
print("Enter a message   : ",end="")
message=int(input())
print("Message sent :     ",message)
get=ReedSolomonSend(message)
print("Message received : ",ReedSolomonReceive(get))
