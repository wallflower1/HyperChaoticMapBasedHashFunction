import sys, collections
import csv

def fy(y0, t, params):

	# ytemp = [0.0 for x in range(len(y0))]
	# ytemp[0] = params[0]*(y0[1]-y0[0])
	# ytemp[1] = y0[0]*(params[2]-y0[2])-y0[1]
	# ytemp[2] = (y0[0]*y0[1])-(params[1]*y0[2])

	ytemp = [0.0 for x in range(len(y0))]
	ytemp[0] = params[0]*(y0[1]-y0[0])
	ytemp[1] = y0[0]*(params[2]-y0[2])-y0[1] + y0[3]
	ytemp[2] = (y0[0]*y0[1])-(params[1]*y0[2])
	ytemp[3] = -params[3]*y0[0]

	return ytemp

def rungeStep(t0, y0, dt, params):
	f = [[0.0 for x in range(len(y0))] for y in range(4)]
	ytemp = [0.0 for x in range(len(y0))]

	ytemp = fy(y0, t0, params)
	f[0] = [dt*i for i in ytemp]
	ynext = [i+j/2.0 for i,j in zip(y0, f[0])]	
	
	ytemp = fy(ynext, t0+dt/2.0, params)
	f[1] = [dt*i for i in ytemp]
	ynext = [i+j/2.0 for i,j in zip(y0, f[1])]

	ytemp =  fy(ynext, t0+dt/2.0, params)
	f[2] = [dt*i for i in ytemp]
	ynext = [i+j for i,j in zip(y0, f[2])]

	ytemp = fy(ynext, t0+dt, params)
	f[3] = [dt*i for i in ytemp]

	y1 = [i+((f0+2.0*f1+2.0*f2+f3)/6) for i,f0,f1,f2,f3 in zip(y0,f[0],f[1],f[2],f[3])]
	return y1

def rungeit(params, y0, t0, N, dt):
	t_out = [0.0 for x in range(N+2)]
	t_out[0] = t0
	y_out = [[0.0 for x in range(len(y0))] for y in range(N+2)]
	y_out[0] = y0
	with open('lorenz.csv', 'w') as f:
		writer = csv.writer(f)
		writer.writerow((y_out[0][0],y_out[0][1],y_out[0][2],y_out[0][3]))
		for x in range(1,N+2):
			#print x
			t_out[x] = t_out[x-1] + dt
			y_out[x] = rungeStep(t_out[x-1], y_out[x-1], dt, params)
			writer.writerow((y_out[x][0],y_out[x][1],y_out[x][2],y_out[x][3]))
	return y_out[N+1]	




def RungaKuttaLorenz(x, y, z, u, k):
	# finding values of x, y, z, u 
	# k : control parameter
	# N : number of iterations
	# a, b, c are parameters set to default 

	params = [10.0, 8.0/3.0, 28.0, k]
	y0 = [x, y, z, u]
	t0 = 0.0
	tn = 30.0
	N = 450
	dt = tn/float(N)

	nth_iter = rungeit(params, y0, t0, N, dt)
	#print nth_iter
	return nth_iter		