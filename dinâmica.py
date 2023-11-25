# --------------------------- Configurações -----------------------------
pi = 3.141592653589793

l1			=	0.8				#[m]
l2			=	0.7				#[m]
m1			=	120.53			#[kg]
m2			=	105.36			#[kg]
t1			=	2*pi			#[rad]
t2			=	2*pi			#[rad]
dt1			=	0				#[rad/s]
dt2			=	0				#[rad/s]
ddt1_var	=	0				#[rad/s^2]
ddt2_var	=	0				#[rad/s^2]
g			=	9.807			#[m/s^2]
I1			=	44.997      	#[kg*m^2]
I2			=	30.115 	    	#[kg*m^2]


tempo_final	=	40	 	#[s]

# Nit é o número de iterações. Não tem unidade física.
Nit			=	1000000	#[]

# Método numérico usado para resolver as equações diferenciais
# Métodos disponíveis: RK4, euler1
metodo				=	"RK4"

# Opção para plotar resultados comparando resultado numérico com experimental
# Caso True, será necessário a inserção de um arquivo com os dados do pêndulo simulado
comparar_com_pendulo_experimental	=	False
arquivo_de_dados					=	"Experimentos_PenDup/PenDup_30_30"
separador							=	"\t\t"

# Seleciona a unidade dos resultados angulares
# Opções: °, rad
unidades 							=	"°"

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
resultados[7].append(l1*sin(t1) + l2*sin(t2))
resultados[8].append(-l1*cos(t1) - l2*cos(t2))
resultados[9].append(0)
resultados[10].append(0)

def paraGrausList(lista):
	return list(map(lambda x: paraGraus(x), lista))

def paraGraus(angulo):
	return ((180/pi) * angulo)


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

def ddt1(t, t1, t2, dt1, dt2, ddt2):
	return -((2 * l1 * l2 * m2 * sin(t1) * sin(t2) + 2 * l1 * l2 * m2 * cos(t1) * cos(t2)) * ddt2 + (2 * l1 * l2 * m2 * sin(t1) * cos(t2) - 2 * l1 * l2 * m2 * cos(t1) * sin(t2)) * dt2**2 + (4 * g * l1 * m2 + 2 * g * l1 * m1) * sin(t1)) / (4 * l1**2 * m2 + l1**2 * m1 + 4 * I1)
	
def ddt2(t, t1, t2, dt1, dt2, ddt1):
	return -((2 * l1 * l2 * m2 * sin(t1) * ddt1 +  2 * l1 * l2 * m2 * cos(t1) * dt1**2 + 2 * g * l2 * m2) * sin(t2) + (2 * l1 * l2 * m2 * cos(t1) * ddt1 - 2 * l1 * l2 * m2 * sin(t1) * dt1**2) * cos(t2)) /(l2**2 * m2 + 4 * I2)

# Aplicação de Runge Kutta de 4ª ordem
def RK4(fun, t, r):

	t1 = r[0]
	t2 = r[1]
	dt1	= r[2]
	dt2 = r[3]
	ddt = r[4]
	yn = r[5]
	
	k1 = passo_de_tempo*fun(t, t1, t2, dt1, dt2, ddt)
	k2 = passo_de_tempo*fun(t + passo_de_tempo/2, t1 + 1/2*k1, t2 + 1/2*k1, dt1 + 1/2*k1, dt2 + 1/2*k1, ddt + 1/2*k1)
	k3 = passo_de_tempo*fun(t + passo_de_tempo/2,  t1 + 1/2*k2, t2 + 1/2*k2, dt1 + 1/2*k2, dt2 + 1/2*k2, ddt + 1/2*k2)
	k4 = passo_de_tempo*fun(t + passo_de_tempo, t1 + k3, t2 + k3, dt1 + k3, dt2 + k3, ddt + k3)
	
	return yn + (1/6)*(k1 + 2*k2 + 2*k3 + k4)


# Resolver por Runge-Kutta de 4ª ordem
def solver_runge_kutta4():
	
	t = 0
	
	global t1, t2, dt1, dt2, ddt1_var, ddt2_var, resultados, arquivo_de_dados, separador, g, energia, t1_travado, t2_travado
	
	# Loop de iteração
	for i in range(Nit):
		t += passo_de_tempo
		
		
		dt1 = RK4(ddt1, t, [resultados[1][i], resultados[2][i], resultados[3][i], resultados[4][i], resultados[10][i], dt1])
		dt2 = RK4(ddt2, t, [resultados[1][i], resultados[2][i], resultados[3][i], resultados[4][i], resultados[9][i], dt2])

		t1 += passo_de_tempo * dt1
		t2 += passo_de_tempo * dt2
		
		resultados[0].append(t)
		resultados[1].append(t1)
		resultados[2].append(t2)
		resultados[3].append(dt1)
		resultados[4].append(dt2)
		resultados[5].append(l1*sin(t1))
		resultados[6].append(-l1*cos(t1))
		resultados[7].append( l1*sin(t1) + l2*sin(t2) )
		resultados[8].append(-l1*cos(t1) - l2*cos(t2))
		resultados[9].append(ddt1_var)
		resultados[10].append(ddt2_var)
	
	arquivo = []
	if comparar_com_pendulo_experimental == True:
		#			 t  x1   y1  x2  y2  t1  t2 dt1 dt2
		#arquivo = [ [], [], [], [], [], [], [], [], []  ]
		
		arquivo = lerSV(arquivo_de_dados, sep=separador)
		
	return resultados, arquivo


# Resolver por Euler (ou runge kutta de 1ª ordem)	
def solver_euler1():
	t=0
	
	res_exp = []
	
	global t1, t2, dt1, dt2, ddt1_var, ddt2_var, resultados, arquivo_de_dados, separador, g
	
	a = g
	
	# Loop de iteração
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
		resultados[7].append( l1*sin(t1) + l2*sin(t2) )
		resultados[8].append(-l1*cos(t1) - l2*cos(t2))
		resultados[9].append(ddt1_var)
		resultados[10].append(ddt2_var)	

	arquivo = []
	if comparar_com_pendulo_experimental == True:
		#			 t  x1   y1  x2  y2  t1  t2 dt1 dt2
		#arquivo = [ [], [], [], [], [], [], [], [], []  ]
		
		arquivo = lerSV(arquivo_de_dados, sep=separador)
		
	return resultados, arquivo


# Plota os gráficos
def graph(res0, res1):
	# Plotando o primeiro pêndulo
	# Estilos úteis: "grayscale", "tableau-colorblind10"
	estilo = "tableau-colorblind10"
	
	global energia
	
	# Define as cores das linhas
	cor_teorico = "blue"
	cor_experimental = "red"
	
	# Configura as unidades de medida angulares:
	if unidades == "°":
		res0[1] = paraGrausList(res0[1])
		res0[2] = paraGrausList(res0[2])
		res0[3] = paraGrausList(res0[3])
		res0[4] = paraGrausList(res0[4])
	
		if comparar_com_pendulo_experimental == True:
			res1[5] = paraGrausList(res1[5])
			res1[6] = paraGrausList(res1[6])
			res1[7] = paraGrausList(res1[7])
			res1[8] = paraGrausList(res1[8])
	
	plt.style.use(estilo)

	plt.figure(figsize=(16,12))

	#	------------------------------------------
	plt.subplot(2, 3, 1)
	plt.plot(res0[0], res0[1], color=cor_teorico)
	plt.grid()
	
	plt.xlabel("t [s]")
	plt.ylabel("θ [" + unidades + "]")
	plt.title("θ_1 em função do tempo")
	
	#	------------------------------------------
	plt.subplot(2,3,2)
	plt.grid()
	
	plt.plot(res0[5], res0[6], color=cor_teorico)

	plt.xlabel("x [m]")
	plt.ylabel("y [m]")
	plt.title("posição da extremidade do primeiro corpo")
	
	#	------------------------------------------
	
	plt.subplot(2,3,3)
	plt.grid()
	plt.plot(res0[1], res0[3], color=cor_teorico)
	
	plt.xlabel("θ [" + unidades + "]")
	plt.ylabel("θ ponto [" + unidades + "/s]")
	plt.title("velocidade angular em função do ângulo do primeiro corpo")
	
	#	------------------------------------------
	plt.subplot(2, 3, 4)
	plt.plot(res0[0], res0[2], color=cor_teorico)
	plt.grid()
	
	plt.xlabel("t [s]")
	plt.ylabel("θ [" + unidades + "]")
	plt.title("θ_2 em função do tempo")
	
	#	------------------------------------------
	plt.subplot(2,3,5)
	plt.grid()
	plt.plot(res0[7], res0[8], color=cor_teorico)
	
	plt.xlabel("x [m]")
	plt.ylabel("y [m]")
	plt.title("posição da extremidade do segundo corpo")
	
	#	------------------------------------------
	
	plt.subplot(2,3,6)
	plt.grid()
	plt.plot(res0[2], res0[4], color=cor_teorico)
	
	plt.xlabel("θ [" + unidades + "]")
	plt.ylabel("θ ponto [" + unidades + "/s]")
	plt.title("velocidade angular em função do ângulo do segundo corpo")
		
	plt.tight_layout()
	plt.show()
	
	# Comparando com o pêndulo experimental, caso dejesado
	if comparar_com_pendulo_experimental == True:
		plt.style.use(estilo)

		plt.figure(figsize=(16,12))

		lim = 0
		# Limita o tempo de comparação do pendulo experimental
		for i in range(len(res1[0])):
			if res1[0][i] >= tempo_final:
				lim = i
				break

		#	------------------------------------------
		plt.subplot(2, 3, 1)
		plt.plot(res1[0][:lim], res1[5][:lim], color=cor_experimental)
		plt.grid()
		
		
		plt.xlabel("t [s]")
		plt.ylabel("θ [" + unidades + "]")
		plt.title("θ_1 em função do tempo")
		
		#	------------------------------------------
		plt.subplot(2,3,2)
		plt.grid()
		plt.plot(res1[1][:lim], res1[2][:lim], color=cor_experimental)
		
		plt.xlabel("x [m]")
		plt.ylabel("y [m]")
		plt.title("posição da extremidade do primeiro corpo")
		
		#	------------------------------------------
		
		plt.subplot(2,3,3)
		plt.grid()
		plt.plot(res1[5][:lim], res1[7][:lim], color=cor_experimental)	
		
		plt.xlabel("θ [rad]")
		plt.ylabel("θ ponto [" + unidades + "/s]")
		plt.title("velocidade angular em função do ângulo do primeiro corpo")
		
		#	------------------------------------------
		plt.subplot(2, 3, 4)
		plt.plot(res1[0][:lim], res1[6][:lim], color=cor_experimental)
		plt.grid()
		
		plt.xlabel("t [s]")
		plt.ylabel("θ [" + unidades + "]")
		plt.title("θ_2 em função do tempo")
		
		#	------------------------------------------
		plt.subplot(2,3,5)
		plt.grid()
		plt.plot(res1[3][:lim], res1[4][:lim], color=cor_experimental)
		
		plt.xlabel("x [m]")
		plt.ylabel("y [m]")
		plt.title("posição da extremidade do segundo corpo")
		
		#	------------------------------------------
		
		plt.subplot(2,3,6)
		plt.grid()
		plt.plot(res1[6][:lim], res1[8][:lim], color=cor_experimental)
		
		plt.xlabel("θ [" + unidades + "]")
		plt.ylabel("θ ponto [" + unidades + "/s]")
		plt.title("velocidade angular em função do ângulo do segundo corpo")
		
		plt.tight_layout()
		plt.show()	
		
		
# Função "Principal"
def main():

	if unidades not in ["°", "rad"]:
		print("Unidade angular inválida")
		return 
	
	if metodo == "RK4":
		res0, res1 = solver_runge_kutta4()
		return graph(res0, res1)
	
	
	if metodo == "euler1":	
		res0, res1 = solver_euler1()
		return graph(res0, res1)
	
	print("Nenhum método válido informado!")
		
main()
