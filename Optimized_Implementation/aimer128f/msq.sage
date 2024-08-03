n = 128
f128 = x^n+x^7+x^2+x+1
F.<X> = GF(2**n, f128)
toHex = lambda x: "{num:0{width}x}".format(num=x.integer_representation(), width=(n/8)*2)
toU64 = lambda x: "0x{:016x}".format(x)
toC = lambda x: ",".join([ toU64(xi) for xi in ZZ(x.integer_representation()).digits(2**64) ])
fromHex = lambda x: F.fetch_int(int(x, 16))
a = F.random_element()
s = toHex(a) # returns a string
s
b = fromHex(s)
a == b # true


def createTable(y, k):
    T = [0]*ceil(n/y)

    for j in range(ceil(n/y)):
        T[j] = [F(0)] * (2**y)
        for i in range(2**y):
            T[j][i] = F(sum([ bit*X**(y*j+z) for (z,bit) in enumerate(ZZ(i).digits(2)) ])) ** (2**k)

    return T

def evalExp(T,y,a):
    c = F(0)
    a_int = a.integer_representation()
    for j in range(ceil(n/y)):
        i = ( a_int >> (y*j) ) % (2**y)
        c += T[j][i]
    return c

k = 91 # exponent
y = 1  # m-ary
T = createTable(y, k)
sizeT = len(T)*len(T[0])
print(sizeT)
a = F.random_element()
b = evalExp(T, y, a)
c = a**(2**k)
print(toHex(a))
print(toHex(b))
print(toHex(c))

def printTableC():
    print( ",\n".join([ "{{{}}}".format( toC(ti[1]) ) for ti in T ]) )
