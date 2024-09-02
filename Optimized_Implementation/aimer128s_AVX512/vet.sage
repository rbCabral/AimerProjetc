F = GF(2)
load("./output.sage")

c_to_matrix = lambda M : matrix([[ F(int(j)) for j in list(i) ] for i in M])
c_to_vet    = lambda M : matrix([[ F(int(j)) for j in list(i) ] for i in M])

vet_b = c_to_vet(VET_A)
vet_c = c_to_vet(VET_C)

mat = c_to_matrix(MAT_B)
vet_R = c_to_vet(VET_R)
vet_R_AVX = c_to_vet(VET_R_AVX)

mat_T = mat.transpose()
AT = vet_c + (vet_b * mat)

print(AT)
print(vet_R_AVX)

assert(AT == vet_R)

assert(vet_R_AVX == vet_R)

