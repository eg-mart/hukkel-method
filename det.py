import numpy as np
import sympy as sp
import pysmiles
import networkx as nx
from sympy.solvers import solve

class polynome:
	def __init__(self, pol):
		self.polynome = pol
	
	def __mul__(self, other):
		polynome_other = ''
		if not self.polynome.startswith('('):
			self.polynome = '(' + self.polynome + ')'
		if not other.polynome.startswith('('):
			polynome_other = '(' + other.polynome + ')'
		else:
			polynome_other = other.polynome
		return polynome(self.polynome + ' * ' + polynome_other)
	
	def __add__(self, other):
		return polynome(self.polynome + ' + ' + other.polynome)
		
	def __sub__(self, other):
		return polynome(self.polynome + ' - ' + other.polynome)


def determinant(matrix):
	result = 0
	if len(matrix) > 3:
		for j in range(1, len(matrix)):
			result += (-1) ** (1 + j) * matrix[1][j] * minor(matrix, 1, j)
		return result
	else:
		return matrix[1][1] * matrix[2][2] - matrix[1][2] * matrix[2][1]


def minor(matrix, i, j):
	new_matrix = [0]
	for m in range(1, len(matrix)):
		if m != i:
			new_matrix.append([0])
			for l in range(1, len(matrix)):
				if l != j:
					new_matrix[-1].append(matrix[m][l])
	return determinant(new_matrix)
	

def generate_ham(bonds, n):
	H = [0]
	a = sp.Symbol('a')
	b = sp.Symbol('b')
	e = sp.Symbol('e')
	
	for i in range(1, n + 1):
		H.append([0])
		for j in range(1, n + 1):
			if (i - 1, j - 1) in bonds or (j - 1, i - 1) in bonds:
				H[-1].append(b)
			elif i == j:
				H[-1].append(a)
			else:
				H[-1].append(0)
	return H


def millennial_equation(H, S):
	e = sp.Symbol('e')
	
	result = [0]
	
	for i in range(1, len(H)):
		result.append([0])
		for j in range(1, len(H)):
			result[-1].append(H[i][j] - S[i][j] * e)
	return result


def generate_system(H, S, E, C):
	equations = []
	for k in range(len(E)):
		cur_equations = []
		for i in range(1, len(H)):
			cur_equation = None
			for j in range(1, len(H)):
				if cur_equation is None:
					cur_equation = C[k + 1][j] * (H[i][j] - E[k])
				else:
					cur_equation += C[k + 1][j] * (H[i][j] - E[k])
			cur_equations.append(cur_equation)
		equations += cur_equations
	last_equation = None
	for i in range(1, len(C)):
		for j in range(1, len(C)):
			if last_equation is None:
				last_equation = C[i][j] ** 2
			else:
				last_equation += C[i][j] ** 2
	last_equation -= 1
	equations.append(last_equation)
	return equations


e = sp.Symbol('e')
smiles = 'C=C'
mol = pysmiles.read_smiles(smiles)
n = len(mol.nodes())

H = generate_ham(mol.edges(), n)
S = [[1 if i == j else 0 for i in range(n + 1)] for j in range(n + 1)]

energy_equation = determinant(millennial_equation(H, S))
E = list(sp.solveset(energy_equation, e))

C_symb = [0]
for i in range(n):
	C_symb.append([0])
	for j in range(n):
		C_symb[-1].append(sp.Symbol('c' + str(i) + str(j)))

C_unknown = []
for i in range(1, len(C_symb)):
	for j in range(1, len(C_symb)):
		C_unknown.append(C_symb[i][j])

C = [0]
solution = solve(generate_system(H, S, E, C_symb), C_unknown)
print(solution)
for i in range(1, n + 1):
	C.append([0])
	for j in range(1, n + 1):
		C[-1].append(solution[C_symb[i][j]])
	print(*C[-1])

print(E)
