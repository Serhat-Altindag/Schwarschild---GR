# -*- coding: utf-8 -*-
"""
Created on Mon Dec 28 00:25:11 2020

@author: MONSTER
"""

import sympy
from einsteinpy.symbolic import MetricTensor, ChristoffelSymbols 
from einsteinpy.symbolic import RiemannCurvatureTensor, RicciTensor, RicciScalar

sympy.init_printing() # for symbolic print out

syms = sympy.symbols("t r theta phi")
G, M, c, a = sympy.symbols("G M c a")

# using metric values of schwarschild space-time
# a is schwarzschild radius; a=2GM/c^2

# 4 tane 1x4luk vektorden 4x4 matrix
matrix4x4 = [[0 for i in range(4)] for i in range(4)]

# Schwarschild metriğin elemanları 
matrix4x4[0][0] = 1 - (a / syms[1])
matrix4x4[1][1] = -1 / ((1 - (a / syms[1])) * (c ** 2))
matrix4x4[2][2] = -1 * (syms[1] ** 2) / (c ** 2)
matrix4x4[3][3] = -1 * (syms[1] ** 2) * (sympy.sin(syms[2]) ** 2) / (c ** 2)

#list ---- > tensor
SCH_Metric = MetricTensor(matrix4x4, syms)
print("SC’s metric ")
SCH_Metric.tensor()

RC_Tensor = RiemannCurvatureTensor.from_metric(SCH_Metric)
print("Riemann Curvature Tensor for SC’s metric")
print("note that “a was Schwarzschild Radius”; \na=2GM/c^2")
RC_Tensor.tensor()

Ric_Tensor = RicciTensor.from_metric(SCH_Metric)
print("Ricci Tensor for SC’s metric ")
Ric_Tensor.tensor()

R = RicciScalar.from_riccitensor(Ric_Tensor)
print("Ricci Scalar for SC’s metric ")
R.expr