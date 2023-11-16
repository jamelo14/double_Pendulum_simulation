# --------------------------- Configurações -----------------------------
pi = 3.141592653589793

l1			=	1		#[m]
l2			=	1		#[m]
m1			=	1		#[kg]
m2			=	1		#[kg]
t1			=	pi/6	#[rad]
t2			=	pi/6	#[rad]
dt1			=	0		#[rad/s]
dt2			=	0		#[rad/s]
ddt1_var	=	0		#[rad/s^2]
ddt2_var	=	0		#[rad/s^2]
g 			=	9.81	#[m/s^2]
I1			=	1		#[kg*m^2]
I2			=	1		#[kg*m^2]

tempo_final	=	40	 	#[s]

# Opção para plotar resultados comparando resultado numérico com experimental
# Caso True, será necessário a inserção de um arquivo com os dados do pêndulo simulado
comparar_com_pendulo_experimental	=	False
arquivo_de_dados					=	"pendulo_experimental"
separador							=	"\t\t"

# Nit é o número de iterações. Não tem unidade física.
Nit			=	1000000	#[]


# -------------------------------------------------------------------------

passo_de_tempo	=	tempo_final/Nit

# Importando bibliotecas
from math import *
import matplotlib.pyplot as plt


# Definindo a matriz dos reultados
#				0	 1    2    3    4    5    6    7    8   9     10
#			    t    t1   t2   dt1  dt2  x1   y1   x2   y2  ddt1  ddt2
resultados = [  [],  [],  [],  [],  [],  [],  [],  [],  [], [],   []  ]

# Carregando os dados iniciais:
resultados[0].append(0)
resultados[1].append(t1)
resultados[2].append(t2)
resultados[3].append(dt1)
resultados[4].append(dt2)
resultados[5].append(l1*sin(t1))
resultados[6].append(-l1*cos(t1))
resultados[7].append(l2*sin(t1) + l2*sin(t2))
resultados[8].append(-l2*cos(t1) - l2*cos(t2))
resultados[9].append(0)
resultados[10].append(0)

def lerSV(nome, sep="\t"):
	
	saida = []
	
	arquivo = open(nome, "r")
	
	for i in arquivo.readline().split(sep):
		saida.append([])
	
	for linha in arquivo:
		lista = linha.split(sep)
		
		for i in range(len(lista)):
			saida[i].append(float(lista[i]))
			
	return saida
		
	
def paraGrausList(lista):
	return list(map(lambda x: paraGraus(x), lista))

def paraGraus(angulo):
	return ((180/pi) * angulo)

def ddt1(t, t1, t2, dt1, dt2, ddt2):
	return -((2*l1*l2*m2*sin(t1)*sin(t2)+2*l1*l2*m2*cos(t1)*cos(t2)+4*I2)*ddt2 + ( 2*l1*l2*m2*sin(t1)*cos(t2) - 2*l1*l2*m2*cos(t1)*sin(t2))*dt2**2 + (4*g*l1*m2+2*g*l1*m1)*sin(t1)) / (4*l1**2 * m2 + l1**2 * m1+4*I2+4*I1)
	
def ddt2(t, t1, t2, dt1, dt2, ddt1):
	return -( (2*l1*l2*m2*sin(t1)*ddt1 + 2*l1*l2*m2*cos(t1)*dt1**2 + 2*g*l2*m2)*sin(t2) + (2*l1*l2*m2*cos(t1)*ddt1 - 2*l1*l2*m2*sin(t1)*dt1**2 ) * cos(t2) + 4*I2*ddt1) / (l2**2 * m2 + 4*I2)

def solver_integracao_numerica():
	t=0
	
	res_exp = []
	
	global t1, t2, dt1, dt2, ddt1_var, ddt2_var, resultados, arquivo_de_dados, separador
	
	for i in range(Nit):
		t += passo_de_tempo
		
		ddt1_var = ddt1(t, t1, t2, dt1, dt2, resultados[10][i-1])
		ddt2_var = ddt2(t, t1, t2, dt1, dt2,  resultados[9][i-1])
		
		dt1 += passo_de_tempo * ddt1_var
		dt2 += passo_de_tempo * ddt2_var
		
		t1 += passo_de_tempo * dt1
		t2 += passo_de_tempo * dt2
		
		resultados[0].append(t)
		resultados[1].append(t1)
		resultados[2].append(t2)
		resultados[3].append(dt1)
		resultados[4].append(dt2)
		resultados[5].append(l1*sin(t1))
		resultados[6].append(-l1*cos(t1))
		resultados[7].append( l2*sin(t1) + l2*sin(t2) )
		resultados[8].append(-l2*cos(t1) - l2*cos(t2))
		resultados[9].append(ddt1_var)
		resultados[10].append(ddt2_var)	


	if comparar_com_pendulo_experimental == True:
		#			 t  x1   y1  x2  y2  dx1 dy1 dx2 dy2
		#arquivo = [ [], [], [], [], [], [], [], [], []  ]
		
		arquivo = lerSV(arquivo_de_dados, sep=separador)
		
		#		    t   t1  t2  dt1  dt2  x1  y1  x2 y2
		res_exp = [ [], [], [], [],  [] 		  		]
		
		for i in range(len(arquivo[0])):
			# adicionando tempo:
			res_exp[0].append(arquivo[0][i])
			
			# adicionando t1 experimental:
			res_exp[1].append(asin(arquivo[1][i]/l1))
		
			# adicionando t2 experimental:
			res_exp[2].append(asin((arquivo[3][i] - arquivo[1][i])/l2))
		
			#adicionando dt1 experimental
			res_exp[3].append(arquivo[5][i]/(l1* cos(res_exp[1][i])))
		
			#adicionando dt2 experimental
			res_exp[4].append( (arquivo[7][i] - arquivo[5][i]) / (l1*cos(res_exp[1][i])) )
			
			#adicionando a posição
			res_exp.append(arquivo[1:5])
		
	graph(resultados, res_exp)

def graph(res0, res1):
	# outros estilos: "grayscale", 
	plt.style.use("tableau-colorblind10")

	plt.figure(figsize=(16,12))

	#	------------------------------------------
	plt.subplot(2, 3, 1)
	plt.plot(res0[0], res0[1])
	plt.grid()
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[0], res1[1])
	
	plt.xlabel("t [s]")
	plt.ylabel("θ [rad]")
	plt.title("θ_1 em função do tempo")
	
	#	------------------------------------------
	plt.subplot(2,3,2)
	plt.grid()
	plt.plot(res0[5], res0[6])
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[5], res1[6])
	
	plt.xlabel("x [m]")
	plt.ylabel("y [m]")
	plt.title("posição da extremidade do primeiro corpo")
	
	#	------------------------------------------
	
	plt.subplot(2,3,3)
	plt.grid()
	plt.plot(res0[1], res0[3])
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[1], res1[3])
	
	plt.xlabel("θ [rad]")
	plt.ylabel("θ ponto [rad/s]")
	plt.title("velocidade angular em função do ângulo do primeiro corpo")
	
	#	------------------------------------------
	plt.subplot(2, 3, 4)
	plt.plot(res0[0], res0[2])
	plt.grid()
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[0], res1[2])
	
	plt.xlabel("t [s]")
	plt.ylabel("θ [rad]")
	plt.title("θ_2 em função do tempo")
	
	#	------------------------------------------
	plt.subplot(2,3,5)
	plt.grid()
	plt.plot(res0[7], res0[8])
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[7], res1[8])
	
	plt.xlabel("x [m]")
	plt.ylabel("y [m]")
	plt.title("posição da extremidade do segundo corpo")
	
	#	------------------------------------------
	
	plt.subplot(2,3,6)
	plt.grid()
	plt.plot(res0[2], res0[4])
	
	if comparar_com_pendulo_experimental == True:
		plt.plot(res1[2], res1[4])
	
	plt.xlabel("θ [rad]")
	plt.ylabel("θ ponto [rad/s]")
	plt.title("velocidade angular em função do ângulo do segundo corpo")
	

	plt.tight_layout()
	plt.show()

solver_integracao_numerica()
