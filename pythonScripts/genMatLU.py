def toHex(h):
    # h = h[::-1]
    b = ''.join(str(bit) for bit in h)
    decimal = int(b, 2)
    hexadecimal = hex(decimal)
    return hexadecimal

def matTrans(v):
    mat = [[0 for x in range(128)] for y in range(128)]
    for i in range(128):
        for j in range(128):
            mat[i][j] = v[j][i]
    return mat

def MaxtoBin(v):
    k=0
    mat0 = [[0 for x in range(128)] for y in range(128)]
    mat1 = [[0 for x in range(128)] for y in range(128)]
    for i in v:
        h=0
        decimal = int(i, 16)  
        binario = bin(decimal)[2:]
        binario = binario.zfill(128)
        for j in binario:
            mat0[k][h] = j
            h = h+1
        k = k+1

    for i in range(128):
        for j in range(128):
            if(mat0[i][j] == '0'):
                mat1[i][j] = 0
            else:
                mat1[i][j] = 1
    
    for i in range(128):
        mat1[i] = mat1[i][::-1]
    return mat1

def matMul(mat1,mat2):
    out = [[0 for x in range(128)] for y in range(128)]
    for i in range(128):
        for j in range(128):
            for k in range(128):
                out[i][j] ^= mat1[i][k] & mat2[k][j]
    return out


def matMulA(mat1,mat2):
    out = []
    for i in range(128):
        # print(mat2[i])
        tmp = [0 for x in range(128)]
        for j in range(128):
            if(mat2[i][j] == 1):
                # print(mat1[i])
                # a = input()    
                for k in range(128):
                    tmp[k] ^= mat1[j][k]
        # print(tmp)
        # a = input()
        out.append(tmp)
    return out



with open('matL', 'r') as arquivo:
    conteudo = arquivo.read()

l = conteudo.splitlines()

with open('matU', 'r') as arquivo:
    conteudo = arquivo.read()

u = conteudo.splitlines()    


l0 = MaxtoBin(l)
u0 = MaxtoBin(u)

u0 = matTrans(u0)

out = matMul(l0,u0)

for i in out:
    print(toHex(i))

# for i in mat1:
#     i = i[::-1]
#     b = ''.join(str(bit) for bit in i)
#     decimal = int(b, 2)
#     hexadecimal = hex(decimal).zfill(16)
    # print(hexadecimal)
    # print(i)
    # print(i[::-1])

# for i in mat1:
    # i = i[::-1]
    # b = ''.join(str(bit) for bit in i)
    # decimal = int(b, 2)
    # hexadecimal = hex(decimal)
    # print(hexadecimal)
    # print(i)
    # print(i[::-1])



