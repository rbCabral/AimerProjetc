v = ['0' for _ in range(128)]

def toHex(h):
    # h = h[::-1]
    b = ''.join(str(bit) for bit in h)
    decimal = int(b, 2)
    hexadecimal = hex(decimal).zfill(8)
    print(hexadecimal)

def mulT(mat0):
    a = '00110101101100000111111011010001000001001011101000101010011110100110010010111100011111111001100100101100010000001000000011100000'
    c = '00000011111111110000011011001001010111111101110001010001010011010000101011101011000111101111101110100011010101011100011111111000'
    a = a[::-1]
    # c = c[::-1]
    s = [int(bit) for bit in a]
    C = [int(bit) for bit in c]
    t = [0 for _ in range(128)]
    for i in range(128):
        if(s[i] == 1):
            # toHex(mat0[i])
            for j in range(128):
                t[j] = t[j] ^ mat0[i][j]
    # print(t[::-1])
    for j in range(128):
        t[j] = t[j] ^ C[j]
    print(t[::-1])

def mul(mat0):
    a = '00110101101100000111111011010001000001001011101000101010011110100110010010111100011111111001100100101100010000001000000011100000'
    c = '00000011111111110000011011001001010111111101110001010001010011010000101011101011000111101111101110100011010101011100011111111000'
    a = a[::-1]
    # c = c[::-1]
    s = [int(bit) for bit in a]
    C = [int(bit) for bit in c]
    t = [0 for _ in range(128)]
    out = []
    # print(s)
    for i in range(128):
        t = 0
        # print(mat1[i])
        for j in range(128):
            t = (s[j] & mat0[i][j]) ^ t 
        out.append(t)
    print(out)
    # for j in range(128):
    #     out[j] = out[j] ^ C[j]
    # print(out)

v[0]   = '0x00000000000000000000000000000001'
v[1]   = '0xfbcf8c1e442c8cf5b87c1159421de6c0'
v[2]   = '0x4d328e5ae8b1bde5687634c0bd8f66a6'
v[3]   = '0x08974511b147a2fe742a6036d93c2057'
v[4]   = '0x4d758c29eeb484f7c8b21bf16608e4db'
v[5]   = '0x82739c986dfbdb200b5c6d5c43980a3c'
v[6]   = '0x34f149a76cf782b00ace7f98da3711b9'
v[7]   = '0xc2ff5007f8a98c83321995ec53ea9914'
v[8]   = '0x097da6d2e8f7686d939b53119c4b7496'
v[9]   = '0x10f77bb9e7748ed35fb6dd3ca90cff95'
v[10]  = '0xc253f8ea7ac0779a55194932141d0937'
v[11]  = '0xdfef8e3e1b142c07b2a4b4591251916b'
v[12]  = '0x931f7bdb443197a114df24dfc33e1f4f'
v[13]  = '0x24128da6bf057bc8bd4cbe8b919dbb07'
v[14]  = '0xb7330162b6115e901be6a922a8d0d7d4'
v[15]  = '0x930f20cea1c668e0b6d9e6635ec916aa'
v[16]  = '0x60488351c7403436ccbb31a458da0423'
v[17]  = '0x9237f55823767eaeef86b4dbc4263e4d'
v[18]  = '0x967e64f599297f3ce2a0e301bed0748a'
v[19]  = '0x58f5315e0e29e3582fde9314f05105e5'
v[20]  = '0xa305f4f11aaa8ad2c9e5b15be18b7596'
v[21]  = '0x31b050cca997ed24a592cb3563071925'
v[22]  = '0x5904c31aaebea1eda55f9e7374b10af1'
v[23]  = '0xea5142776b77d368cf6921d88d12bbf2'
v[24]  = '0x448bfd74cc62450628779ef24c9ddcb5'
v[25]  = '0xc66ef14828e98e800d2caf1924759d9a'
v[26]  = '0x5121956dac40960a312a49ac8d3790db'
v[27]  = '0x41cdda4642d1e45a311230a0f0166f37'
v[28]  = '0x50accd8f44cc6a3b152cd68f8d980779'
v[29]  = '0xaf2e59b6e13fec010e6342e6e178a202'
v[30]  = '0x008f7a68483d8f8acdfaea274cfff823'
v[31]  = '0xdece92499f9521ba80183f4571309485'
v[32]  = '0x3c5814a4c792b3beaba321469362905f'
v[33]  = '0x7585a167f0b843b87e8680766e1d3ffc'
v[34]  = '0x114bd1d466ef87874e81e572c5dbf79a'
v[35]  = '0x014747267c0b38f93a7e0a403a1da600'
v[36]  = '0x196284a6305e23d723116c4dd539e293'
v[37]  = '0xae9980fc3aa3cd2ef0a02f00d5a45c0f'
v[38]  = '0xd2be72a1da8addbc7eae2c6dae8286e9'
v[39]  = '0x2d1eb9e86163e7f4bb8689cb630a9e23'
v[40]  = '0x8eda5ccd665a3ac2f0febfb8f6e46561'
v[41]  = '0x83fe14fe68c33df0370a6880719f8be9'
v[42]  = '0xdfb51a0ca76c9c82e9634dd58474116f'
v[43]  = '0xcecdf64c8559eef49c40da32ca69fe52'
v[44]  = '0x9256190cf3cfb1fae29f358edce8d40e'
v[45]  = '0xc69025ae5a99210ab5431f672597e9cd'
v[46]  = '0xdf81e3ce7617b0a10f00e0c670d40d95'
v[47]  = '0xc24d5671c235f28e699332d099ea38d7'
v[48]  = '0x7c11f6654369b65d89ea2f4529a74b45'
v[49]  = '0x049b62170967135beaa8470e44915e89'
v[50]  = '0xba3743d76fda508339fb9877aadc951b'
v[51]  = '0xdef2a5ff6e028abae2da8722532e6fbb'
v[52]  = '0xf28418e25fbc0144b5e975340c6c76a1'
v[53]  = '0xab56f227d4a26a26035ab9363f6882be'
v[54]  = '0x75af981a11d43e64273536b8b02dd5f1'
v[55]  = '0x507a048207335fa1846e480a8bc44fa9'
v[56]  = '0xf632f1c9c802ab763808d8fa4fcba922'
v[57]  = '0x2dcfedd3c12f73be34ecb7872eda1962'
v[58]  = '0xf77d23a1c6600553f884a540c1b411ea'
v[59]  = '0x7d5ef83763344eed0e106a0239843e3c'
v[60]  = '0x5070be659c9249dd4192e743be4ae7d9'
v[61]  = '0x09d7b6469e9538566588c07b62dd03ba'
v[62]  = '0x5c859acd40414177790b4af55db42c92'
v[63]  = '0xd728f7e92e3e7940edda860c739ca8f7'
v[64]  = '0xf6fd92c58b52c44dfbcf513b18b860f7'
v[65]  = '0x04286d00eb3471974f1571762119854e'
v[66]  = '0xf68288c09c8d73d43f777b9977ed2aa6'
v[67]  = '0x86437c4cb491c94b538b16a3bd887a20'
v[68]  = '0x97db137363bf2a7e3656d64f9fdf8baf'
v[69]  = '0x213b4a759760ffe00582fbdad31a1e6c'
v[70]  = '0x05cb6fb77aad0666c7f42208feff0a47'
v[71]  = '0xd3740dabc91a5ecd8f59c644fd5259d4'
v[72]  = '0xa2486f3cdc03c63fca19d9ef4ab67cee'
v[73]  = '0x71b6e4a0b3d4a2a08a1f14a7c3d2f88f'
v[74]  = '0xa28d2266c47e02b2e9ee9aa288652690'
v[75]  = '0x689aa81596670031759c7eee1a3eead3'
v[76]  = '0x206b34f2ed6fc8ff50a9a3f15e0032ae'
v[77]  = '0xf7f5952347d531a6630774b85c40302e'
v[78]  = '0x7072ec9b3a2ddd8f78886ef4e794267b'
v[79]  = '0xb360d5ec03ebf053754c7bf46deec1a2'
v[80]  = '0xbe8541809bccdc7d337080ceace4b67a'
v[81]  = '0x87fc6f3fbe554f618c243c5d486009a0'
v[82]  = '0xc7a500e89b1b40a658e8f3ccf2596f26'
v[83]  = '0xe719cb9a5a49f8ed516a6cbee9e76420'
v[84]  = '0x484947f2b48d788296f150816f90c216'
v[85]  = '0x58823337d6c0a641defb92978dfa0053'
v[86]  = '0xdad5891c70205c9598bbc22dd2d3262b'
v[87]  = '0xea3e0a86c4241c1ebf1d06e5edc7d9ba'
v[88]  = '0x1110bc39a35669cf78e2cf480abc18ef'
v[89]  = '0x8eb4cf8cb0851480c188299c1375e7b2'
v[90]  = '0xff5c57071581e3b1d0ec275048c667d2'
v[91]  = '0xcf79008ac79991d3955c8d54a50fdd52'
v[92]  = '0x1c5fc0acfab2cbb2f46cdcd85b7289c9'
v[93]  = '0x862183d1a9042f4d676f48ac3ed3c825'
v[94]  = '0xa655183af862baaef35fc7982c7daee6'
v[95]  = '0x963ed04a2a0b3eaa5335bbcaf8b9f37b'
v[96]  = '0x82f1d3e8253374ee76d009714121cb10'
v[97]  = '0xee023bd013e359f550198339f3198270'
v[98]  = '0x5c1520117e098dbf315d27ea94c7941a'
v[99]  = '0x7d84bbe2e504ff9496ccc513ba987df7'
v[100] = '0x7d9fc4a633f228f403464584b630d2b7'
v[101] = '0x45a5498048f1a4747e39cbb756cac943'
v[102] = '0x4883787b94c9042556a90669f7aa29c6'
v[103] = '0x6495beb53f9054019a262b27cb8de6e9'
v[104] = '0x6c2c9c31b3faab04dc5866e0159b2920'
v[105] = '0x1926807fb1c2bdd582f93c693fec7b5f'
v[106] = '0xff56ec036c5f13d63a06ca560fda4251'
v[107] = '0xaea98fd960fd6b9dcf96fe4ae095a1c3'
v[108] = '0xe7c1f21b63d4e19cc2ae3b23e1b73447'
v[109] = '0xa61e4689ac8893c0660f92196e62044c'
v[110] = '0xb71adc881811c2584aacc983cc5d9cfe'
v[111] = '0x3d4b38fea83810f8b01938e5f92ea2e0'
v[112] = '0x242e99777aeec42a8195527abb10f039'
v[113] = '0x928620c22f148a6d077a36f6536baf7f'
v[114] = '0x300ecf50c00b75a7b4d16665e8f965a6'
v[115] = '0x50967628985e6eaf53d4fbf144350d5c'
v[116] = '0x74a182255aa9ccaeea67291009feb48e'
v[117] = '0xbe1b4991d245fa61e67c52e63c97fb3a'
v[118] = '0xc26bdd871e8691e56bd8d3685ed38551'
v[119] = '0x2408058c7b3e3c09267c4e3df39e0a7e'
v[120] = '0xe70baa2724b374d32bc55550057b4b4e'
v[121] = '0x4f4e64ba26bfee680e2984947284c4dd'
v[122] = '0x357f7d8801646f0878891ea4bacdb828'
v[123] = '0x8e6c9653552802fa220a9cb569d1ee6b'
v[124] = '0xcbd0c89374b1cc2e6159359f74dda4d7'
v[125] = '0xeeca37f94d3f69bd8dd5a4c4fe55c89f'
v[126] = '0x69b86caf61d48d1522abf1f68e0f314a'
v[127] = '0x1a49957d9798f177ab26c59f1090455d'

mat0 = [[0 for x in range(128)] for y in range(128)]
mat1 = [[0 for x in range(128)] for y in range(128)]
mat2 = [[0 for x in range(128)] for y in range(128)]

k=0
#ler a matriz
for i in v:
    h=0
    decimal = int(i, 16)  
    binario = bin(decimal)[2:]
    # print(binario)
    binario = binario.zfill(128)
    for j in binario:
        mat0[k][h] = j
        h = h+1
    k = k+1

#passa para inteiro
for i in range(128):
    for j in range(128):
        if(mat0[i][j] == '0'):
            mat0[i][j] = 0
        else:
            mat0[i][j] = 1

mulT(mat0)


#retornar transposta
for i in range(128):
    for j in range(128):
        mat1[i][j] = mat0[j][i]

for i in range(128):
    mat1[i][0:64]     = mat1[i][63::-1]
    mat1[i][64:128]   = mat1[i][127:63:-1]

for i in range(16):
    for j in range(8):
        mat2[i*8+j] = mat1[i*8+7-j]
        # print(i*8+j,i*8+7-j)


# mul(mat1)

# for i in mat2:
#     toHex(i)

# a = '00110101101100000111111011010001000001001011101000101010011110100110010010111100011111111001100100101100010000001000000011100000'
# c = '00000011111111110000011011001001010111111101110001010001010011010000101011101011000111101111101110100011010101011100011111111000'
# a = a[::-1]
# c = c[::-1]
# s = [int(bit) for bit in a]
# C = [int(bit) for bit in c]
# # t = [0 for _ in range(128)]
# out = []
# for i in range(128):
#     t = 0
#     for j in range(128):
#         t = (s[j] & mat1[i][j]) ^ t 
#     out.append(t)

# print(out)

# # for i in mat1:
# #     print(i[::-1])



# # for i in mat1:
# #     i = i[::-1]
# #     b = ''.join(str(bit) for bit in i)
# #     decimal = int(b, 2)
# #     hexadecimal = hex(decimal).zfill(8)
# #     print(hexadecimal)
#     # print(i)
#     # print(i[::-1])

# # for i in mat1:
#     # i = i[::-1]
#     # b = ''.join(str(bit) for bit in i)
#     # decimal = int(b, 2)
#     # hexadecimal = hex(decimal)
#     # print(hexadecimal)
#     # print(i)
#     # print(i[::-1])



