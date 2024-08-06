F = GF(2)
load("./output.sage")

c_to_matrix = lambda M : matrix([[ F(int(j)) for j in list(i) ] for i in M])

c_LT = c_to_matrix(LT)
c_UT = c_to_matrix(UT)
c_AT = c_to_matrix(AT)

AT = c_UT * c_LT

assert(AT == c_AT)  # true

_,ll,uu = AT.LU()
assert(ll == c_UT)
assert(uu == c_LT)

# Verification
A = AT.transpose()
L = c_LT.transpose()
U = c_UT.transpose()
assert( A == L*U ) 


