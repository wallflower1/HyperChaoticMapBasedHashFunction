import collections
import sys, os, binascii
import math, random, operator
import RK4lorenz

# step1: take msg input, parse it into binary representation
# step2: pad the msg so that its length is a multiple of M(128)
# step3: generate key parameters h1,h2,h3,h* of length N(128)
# step4: iterate the initial values through the pipe structure, storing intermediate values
# step5: generate the final hash value (128 bit) and store in a file

def compression_fn(p1, p2, p3):
	p = [[0 for x in range(4)] for y in range (3)]

	for y in xrange(0,4):
		p[0][y] = int(p1[y*32:y*32+31],2)
	for y in xrange(0,4):
		p[1][y] = int(p2[y*32:y*32+31],2)
	for y in xrange(0,4):
		p[2][y] = int(p3[y*32:y*32+31],2)
	
	x0 = (operator.xor(p[2][0], p[1][0]) + operator.xor(p[0][0], p[2][0]))/math.pow(2,32)
	y0 = (operator.xor(p[2][1], p[1][1]) + operator.xor(p[0][1], p[2][1]))/math.pow(2,32)
	z0 = (operator.xor(p[2][2], p[1][2]) + operator.xor(p[0][2], p[2][2]))/math.pow(2,32)
	u0 = (operator.xor(p[2][3], p[1][3]) + operator.xor(p[0][3], p[2][3]))/math.pow(2,32)
	k = (operator.xor(p[2][1],p[1][1]) + operator.xor(p[0][1],p[2][1]))/math.pow(2,28)

	if k<0 or k>16.6:
	 	k = 5
	#print x0,y0,z0,u0,k

	nth_iter = RK4lorenz.RungaKuttaLorenz(x0, y0, z0, u0, k)
	for x in xrange(0,4):
		val = abs(nth_iter[x])
		while val < math.pow(2,31)-1:
			val *= 10
		val /= 10
		nth_iter[x] = str(bin(int(math.floor(val))))[2:].rjust(32,'0')
	# print nth_iter 
	# print len(nth_iter[0]),len(nth_iter[1]),len(nth_iter[2]),len(nth_iter[3])
	
	h = [[0 for x in range(4)] for y in range (4)]
	for x in xrange(0,4):
		for y in xrange(0,4):
			h[x][y] = nth_iter[x][y*8:y*8+8]	
	#print h
	Q1 = h[0][0]+h[1][0]+h[2][0]+h[3][0]+h[0][2]+h[1][2]+h[2][2]+h[3][2]
	Q2 = h[0][1]+h[1][1]+h[2][1]+h[3][1]+h[0][3]+h[1][3]+h[2][3]+h[3][3]
	Q = Q1+Q2
	#print hex(int(Q,2))
	return Q	


def generate_hash(M, N, msg_blocks, H):
	L = len(msg_blocks)
	P = []
	O = [0 for x in range(L)]
	H_i = H

	for x in xrange(0,L):
		#print H_i
		h1 = compression_fn(H_i[0], H_i[1],msg_blocks[x])
		h2 = compression_fn(H_i[1], H_i[2],msg_blocks[x])
		h3 = compression_fn(H_i[2], H_i[0],msg_blocks[x])
		h4 = compression_fn(H_i[0], H_i[1], H_i[2])
		O[x] = bin(int(h1,2)^int(h2,2)^int(h3,2)^int(h4,2))[2:].rjust(N,'0')
		#print O[x]
		H_i[0] = h1
		H_i[1] = h2
		H_i[2] = h3
	final_O = 0
	for x in xrange(0,L):
		final_O = final_O ^ int(O[x])
	#	print final_O
	final_O = bin(final_O)[2:].rjust(N,'0')

	h1 = compression_fn(H_i[0], H_i[1],final_O)
	h2 = compression_fn(H_i[1], H_i[2],final_O)
	h3 = compression_fn(H_i[2], H_i[0],final_O)
	h4 = compression_fn(H_i[0], H_i[1], H_i[2])

	result_hash = hex(int(h1,2)^int(h2,2)^int(h3,2)^int(h4,2))[2:].rjust(32,'0')
	return result_hash
def generate_params(n):
	
	f = open('init_parameter.txt','w+')
	h = []
	for x in xrange(1,4):
		h.append(bin(random.randint(0,math.pow(2,n))))
	for x in xrange(0,3):
		f.write(str(h[x][2:]).ljust(128,'0')+' ')
	f.close()
	compression_fn(h[0], h[1], h[2])
	return h

def break_msg(msg, N):
	msg_blocks = []
	length = len(msg)
	k=0
	
	while k < length:
		msg_blocks.append(msg[k:k+N])
		k = k+N

	return msg_blocks

# padding the msg block according to Merkel-Damgard padding scheme
def pad_msg(msg, block_size):
	msg_len = len(msg)
	if msg_len % block_size == 0:
		final_msg = msg + '1'.ljust(block_size, '0')
	elif msg_len % block_size < 127:
		msg = msg+'1'
		final_msg = msg + ''.ljust(block_size-(msg_len%block_size)-1, '0')
	else:
		final_msg = msg + '1' + ''.ljust(block_size, '0')
	#print 'padded: '+ final_msg
	#print len(final_msg)
	return final_msg


def construct_hash(msg, hash_file):
	# M: length of individual msg block (in bytes)
	# N: length of final and intermediate hash values (in bytes)
	M = 16 
	N = 16

	#read msg from file
	#message = ''
	#input_msg_file = open('input_msg.txt', 'rb')
	# for x in input_msg_file:
	# 	message = message + x
	msg_bin = str(bin(int(binascii.b2a_hex(msg), 16)))
	
	#pad message
	padded_msg = pad_msg(msg_bin, M*8) 
	
	msg_blocks = break_msg('00'+padded_msg[2:], N*8)
	
	H = [0 for x in range(3)]
	H[0]='10110101100111111001010101111010111110101000111010010010110000111110010101000011001010111110101000011110010000110101010101000011'
	H[1]='10000111101001110000011000001110100000111111111111101011101011101001100101101110101110111110101100011110001101010111001100000100'
	H[2]='10111011010010100010100111001001100001100111001011111011001011100110011110111001110011101100000111101010000101101010000011010011'
	
	final_hash_value = generate_hash(M, N*8, msg_blocks, H)[:-1].ljust(32,'0')
	print 'Hash : '+final_hash_value
	hash_file.write(final_hash_value+'\n')
	print len(final_hash_value)

def performance_check():
	msg = list('this is sample text.')
	with open('hash.txt','w') as hash_file:
		for x in xrange(0,100000):
			y = random.randint(0, len(msg)-1)
			z = random.randint(0, 256)
			msg[y] = chr(97+(ord(msg[y])+z)%26)
			print 'msg: '+''.join(msg)
			construct_hash(''.join(msg), hash_file)


if __name__ == '__main__':
	performance_check()