import mpmath as mp
from flint import arb, acb, acb_poly, ctx, fmpq, good

import psutil
from math import floor
import time
import sys

"""
This python sheet calculates the matrix elements for the C and P matrices, associated with a multiple backflow
operator, as detailed in "Repeated Quantum backflow and overflow" (https://arxiv.org/abs/2505.13184). 

The parameters below mean as follows:
N - Largest trial vector used in the sense that the largest backflow and overflow spectral 
points are found in the subspace of L^2 functions spanned by (psi__n)__(0<=n=N).
M - Number of backflow periods under investigation
delta - Shift of the exponent of the trial vectors as defined above
A - Exponential multiplier of the trial vectors as defined above

***By default the Python code works with the assumption that the backflow times under investigation are given by 
t__k = (2*k - 2*M + 1)/2 for k = 0 .. 2*M - 1. Particular times can be specified using the 
times variable.***
dps - Digits of the file, by default this is set to N.
"""


def check_times(array):
    if len(array) % 2 != 0:
        raise Exception(
            "Times vector must contain an even number of elements.")
    for ii in range(0, len(array)-1):
        if array[ii] > array[ii+1]:
            raise Exception("Times must be in increasing order.")


def to_fmpq(x):
    n = len(str(x))
    return fmpq(int(10**n*float(x)), 10**n)


custom_times = ',' in sys.argv[4]

if custom_times:
    mode = 0
    try:
        N_start = int(sys.argv[1])
    except:
        raise Exception("First argument N_start must be an integer.")
    try:
        N_stop = int(sys.argv[2])
    except:
        raise Exception("Second argument N_stop must be an integer.")
    try:
        N_max = int(sys.argv[3])
    except:
        raise Exception("Third argument N_max must be an integer.")
    try:
        print(sys.argv[4])
        times = [to_fmpq(element) for element in sys.argv[4].split(",")]
        str_times = sys.argv[4]
        M = floor(len(times)/2)
    except:
        raise Exception(
            "Second argument must be a comma separated list of increasing numbers or integer number of backflow periods.")
    if len(times) % 2 != 0:
        raise Exception("Provide an even number of times.")
    check_times(times)
    try:
        str_delta = str(sys.argv[5])
        delta = to_fmpq(str_delta)
    except:
        raise Exception("Third argument delta must be a real number.")
    try:
        str_a = str(sys.argv[6])
        a = to_fmpq(str_a)
    except:
        raise Exception("Fourth argument a must be a positive real number.")
    if a <= 0:
        raise Exception("Fourth argument a must be a positive real number.")
    dps = int(floor(N_max))
    print("Arguments accepted.")
    print("N: ", N_start, "__", N_stop)
    print("Times: ", times)
    print("delta", str_delta)
    print("a", str_a)

elif not custom_times:
    mode = 1
    try:
        N_start = int(sys.argv[1])
    except:
        raise Exception("First argument N_start must be an integer.")
    try:
        N_stop = int(sys.argv[2])
    except:
        raise Exception("Second argument N_stop must be an integer.")
    try:
        N_max = int(sys.argv[3])
    except:
        raise Exception("Third argument N_max must be an integer.")
    try:
        M = int(sys.argv[4])
    except:
        raise Exception(
            "Fourth argument must be a comma separated list of increasing numbers or integer number of backflow periods.")
    try:
        str_delta = str(sys.argv[5])
        delta = to_fmpq(str_delta)
        print(delta)
    except:
        raise Exception("Fifth argument delta must be a real number.")
    try:
        str_a = str(sys.argv[6])
        a = to_fmpq(str_a)
    except:
        raise Exception("Sixth argument a must be a positive real number.")
    if float(str_a) <= 0:
        raise Exception("Sixth argument a must be a positive real number.")
    times = [fmpq(2*ii-2*M+1, 2) for ii in range(0, 2*M)]
    times=[mp.power(mp.pi,-1)*7*M*mp.mpf(-1-2*(M-1)+2*k) for k in range(0,2*M)]
    dps = int(floor(N_max))
    print(times)
    print("Arguments accepted.")
    print("N: ", N_start, "_", N_stop)
    print("M: ", M)
    print("delta", str_delta)
    print("a", str_a)
else:
    print("Command line arguments should not contain more than 6 entries.")
    print("The two acceptable forms are:")
    print("N [t1,t2,...,t2M-1,t2M] delta a")
    print("or")
    print("N M delta a")
    raise Exception(
        "See above for the accepted format of the command line arguments. Please try again.")

ctx.dps = dps+50
mp.mp.dps = dps
print("dps", dps)

I = acb(0, 1)
q = fmpq(1, 4)
h = fmpq(1, 2)

print("a,", a)
print("delta,", delta)

# delta=fmpq(-1,4)
# a=fmpq(636,1000)
# s1,s2,s3,s4=fmpq(-3,2),fmpq(-1,2),fmpq(1,2),fmpq(3,2)
# times=[s1,s2,s3,s4]


def signum(x):
    if x > 0:
        return 1
    elif x < 0:
        return -1
    else:
        return 0

def D(pm, n, m, delta, a):
    mu, nu = acb(n+delta+h), acb(m+delta+h)
    coeff = good(lambda: acb(2*a).pow(n+m+2*delta+1)*acb(n+delta+pm*q+1).gamma()
                 * acb(m+delta-pm*q+1).gamma()*(acb(2*n+2*delta+1).gamma()
                    * acb(2*m+2*delta+1).gamma()).pow(-fmpq(1,2))
                 )
    return coeff


def J(s1, s2, alpha, beta, a):
    if s1*s2 > 0:
        Jval = good(lambda: I*acb(2*a).pow(-alpha-beta-1)*(
            acb(-signum(s1)*(beta+1)).exp_pi_i()
            * (2*a*acb(a-I*s1).pow(-1)).beta_lower(alpha+beta+1, -beta)
            - acb(-signum(s2)*(beta+1)).exp_pi_i()
            * (2*a*(a-I*s2).pow(-1)).beta_lower(alpha+beta+1, -beta)
        )
        )
    elif s1 == -s2:
        Jval = good(lambda: 2*acb(2*a).pow(-alpha-beta-1)*(
            (acb(-signum(s2)*(beta+1)).exp_pi_i()
                * (2*a*(a-I*s2).pow(-1)).beta_lower(alpha+beta+1, -beta)).imag
            + arb.pi()*((alpha+beta+1)*acb(1).beta_lower(alpha+1, beta+1))
            .pow(-1)
        )
        )
    else:
        Jval = good(lambda: I*acb(2*a).pow(-alpha-beta-1)*(
            acb(-signum(s1)*(beta+1)).exp_pi_i()
            * (2*a*acb(a-I*s1).pow(-1)).beta_lower(alpha+beta+1, -beta)
            - acb(-signum(s2)*(beta+1)).exp_pi_i()
            * (2*a*acb(a-I*s2).pow(-1)).beta_lower(alpha+beta+1, -beta)
            - 2*arb.pi()*I *
            ((alpha+beta+1)*(1).beta_lower(alpha+1, beta+1)).pow(-1)
            )
        )
    return Jval
"""

def J(s1, s2, alpha, beta, a):
    if s1*s2 > 0:
        Jval = good(lambda: I*acb(2*a).pow(-alpha-beta-1)*(
            acb(-signum(s1)*(beta+1)).exp_pi_i()
            * acb(2*a).pow(alpha+beta+1)*acb(a,-s1).pow(-alpha-beta-1)
            * acb(alpha+beta+1).pow(-1)
            * acb(2*a*acb(a,-s1).pow(-1)).hypgeom_2f1(alpha+beta+1,1+beta,alpha+beta+2)
            - acb(-signum(s2)*(beta+1)).exp_pi_i()
            * acb(2*a).pow(alpha+beta+1)*acb(a,-s2).pow(-alpha-beta-1)
            * acb(alpha+beta+1).pow(-1)
            * acb(2*a*acb(a,-s2).pow(-1)).hypgeom_2f1(alpha+beta+1,1+beta,alpha+beta+2)
        )
        )
    elif s1 == -s2:
        Jval = good(lambda: 2*acb(2*a).pow(-alpha-beta-1)*(
            (acb(-signum(s2)*(beta+1)).exp_pi_i()
             * acb(2*a).pow(alpha+beta+1)*acb(a,-s2).pow(-alpha-beta-1)
             * acb(alpha+beta+1).pow(-1)
             * acb(2*a*acb(a,-s1).pow(-1)).hypgeom_2f1(alpha+beta+1,1+beta,alpha+beta+2)).imag
             + arb.pi()*((alpha+beta+1)*acb(1).beta_lower(alpha+1, beta+1))
            .pow(-1)
        )
        )
    else:
        Jval = good(lambda: I*acb(2*a).pow(-alpha-beta-1)*(
            acb(-signum(s1)*(beta+1)).exp_pi_i()
            * acb(2*a).pow(alpha+beta+1)*acb(a,-s1).pow(-alpha-beta-1)
            * acb(alpha+beta+1).pow(-1)
            * acb(2*a*acb(a,-s1).pow(-1)).hypgeom_2f1(alpha+beta+1,1+beta,alpha+beta+2)
            - acb(-signum(s2)*(beta+1)).exp_pi_i()
            * acb(2*a).pow(alpha+beta+1)*acb(a,-s1).pow(-alpha-beta-1)
            * acb(alpha+beta+1).pow(-1)
            * acb(2*a*acb(a,-s1).pow(-1)).hypgeom_2f1(alpha+beta+1,1+beta,alpha+beta+2)
            - 2*arb.pi()*I *
            ((alpha+beta+1)*(1).beta_lower(alpha+1, beta+1)).pow(-1)
            )
        )
    return Jval
"""


def matrix_element(n, m, delta, a):
    def D(n, m, delta, a):
        mu, nu = acb(n+delta+h), acb(m+delta+h)
        """
        coeff = good(lambda: acb(2*a).pow(n+m+2*delta+1)*acb(n+delta+pm*q+1).gamma()
                     * acb(m+delta-pm*q+1).gamma()*(acb(2*n+2*delta+1).gamma()
                        * acb(2*m+2*delta+1).gamma()).pow(-fmpq(1,2))
                     )
        """
        coeff=good(lambda: arb(2).sqrt()*arb.pi()
                   *acb(2*a).pow(mu+nu)*(acb(1).beta_lower(mu,mu)
                   *acb(1).beta_lower(nu,nu)).sqrt()
                   *(acb(1).beta_lower(mu,fmpq(3,4))
                     *acb(1).beta_lower(nu,fmpq(1,4))).pow(-1))
        return coeff

    c1, c2 = acb(0),acb(0)
    ctx.dps += 2*M
    for ii in range(0, floor(len(times)/2)):
        c1 += J(times[2*ii], times[2*ii+1],n+delta+q,m+delta-q,a)
        c2 += J(times[2*ii], times[2*ii+1],n+delta-q,m+delta+q,a)
    #p1, p2 = acb_poly(c1), acb_poly(c2)
    # print(D(n,m,delta,a)*p1.evaluate([1])[0])
    # print()
    # print(D(m,n,delta,a)*p2.evaluate([1])[0])
    # print()
    ctx.dps -= 2*M
    return (-(D(n, m, delta, a)*c1
              + D(m, n, delta, a)*c2)
            * acb(4*arb.pi()).pow(-1))


def gram_element(n, m, delta, a):
    mu, nu = acb(n+delta+h), acb(m+delta+h)
    return (good(lambda: (acb(1).beta_lower(mu, mu) *
                 acb(1).beta_lower(nu, nu)).sqrt()
                 * acb(1).beta_lower(mu, nu).pow(-1)
                 )
            )


def convert_to_maple_string(x):
    z = mp.mpmathify(x.str(dps+5, radius=False))
    if mp.im(z) < 0:
        return str(mp.re(z))+str(mp.im(z))+"*I"
    elif mp.im(z) > 0:
        return str(mp.im(z))+"+"+str(mp.im(z))+"*I"
    else:
        return str(mp.re(z))


def data_dict():
    data_dict = dict()
    for n in range(N_start, N_stop+1):
        for m in range(0, n+1):
            if n % 5 == 0 and m % 5 == 0:
                print("dps", mp.mp.dps)
                print("cpu usage: ", psutil.cpu_percent(), "%")
                print("virtual memory available: ",
                      psutil.virtual_memory().available*100 /
                      psutil.virtual_memory().total, "%")
                print('('+str(n)+','+str(m)+')')

            data_dict['matrix_element', n, m] = convert_to_maple_string(
                matrix_element(n, m, delta, a))
            data_dict['gram_element', n, m] = convert_to_maple_string(
                gram_element(n, m, delta, a))

    return data_dict

print("Calculating...")
t0 = time.time()

strindex = ["gram_element", "matrix_element"]
maple_dict = data_dict()

if mode == 0:
    print(str_times)
    out_file = ('auxvals/arb/MB_N'+str(N_start)+'_'+str(N_stop)+'MAX'+str(N_max)
                + '_t'+str_times+'_delta'+str_delta+'_a'+str_a+'_DPS' + str(dps) + '.mpl')
    print(out_file)
elif mode == 1:
    out_file = ('auxvals/arb/MB_N'+str(N_start)+'_'+str(N_stop)+'MAX'+str(N_max)
                + '_M'+str(M)+'_delta'+str_delta+'_a'+str_a+'_DPS' + str(dps) + '.mpl')
else:
    out_file = 'unknown.mpl'


def write_maple(d):
    outlist = []
    for (k, v) in d.items():
        # k is a tuple; its string representation is exactly what is needed
        # as a Maple key
        km = str(k)
        vm = str(v)
        outlist.append(km + ' = ' + vm)
    outstring = 'table([\n' + ',\n'.join(outlist) + '\n])'
    return(outstring)


maple_dict_list = list(maple_dict.items())
maple_dict_list.sort(key=lambda y: y[1], reverse=True)

with open(out_file, 'w') as fo:
    fo.write('data_python := ')
    fo.write(write_maple(maple_dict))
    fo.write(';\n')
print('...finished and saved. Calculation took ', time.time()-t0, 's.')