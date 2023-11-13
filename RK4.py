Nit = 1000
tfinal = 10

h = tfinal / Nit

pi = 3.141592653589793

theta1_ini		=	pi/4
theta2_ini		=	pi/2
dtheta1_ini		=	0
dtheta2_ini		=	0
ddtheta1_ini	=	0
ddtheta2_ini	=	0

#				t	t1	t2	dt1	dt2 ddt1 ddt2
resultados = [  [], [], [], [], []            ]

r		   = [      [], [], [], [], [],  []   ]

resultados[0][0] = 0
resultados[1][0] = theta1_ini
resultados[2][0] = theta2_ini
resultados[3][0] = dtheta1_ini
resultados[4][0] = dtheta2_ini

r[0] = theta1_ini
r[1] = theta2_ini
r[2] = dtheta1_ini
r[3] = dtheta2_ini
r[4] = ddtheta1_ini
r[5] = ddtheta2_ini

def RK4(f, r, t, h):

	k1 = h*f(r,t)
    k2 = h*f(r+k1/2,t+h/2)
    k3 = h*f(r+k2/2,t+h/2)
    k4 = h*f(r+k3,t+h)
    
    return (1/6)*(k1+2*k2+2*k3+k4)
    
	
def f(r, t):
	theta		=	r[0]
	dtheta		=	r[1]
	theta0		=	r[2]
	dtheta0		=	r[3]
	ddtheta0	=	r[4]
	
	return
		 
def sim():
	for t in range(0, tfinal, h):
		theta1 = theta1_ini
		theta2 = theta2_ini
		dtheta1 = 
		
		dtheta1 += RK4(ddt1, r, t, h)
		dtheta2 += RK4(ddt2, r, t, h)
		
		
		
		r[0] = theta1
		r[1] = theta2
		r[2] = dtheta1
		r[3] = dtheta2
		r[4] = ddt1(r, t)
		r[5] = ddt2(r, t)

		
