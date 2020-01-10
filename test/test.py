import toolbox
import numpy as np


#对高斯消去法的测试
print("\n")
print("Test of Gauss Elimination algorithm.")
matrix = np.array([[1, 1, 1, 3],
                   [1, 2, 4, 7],
                   [1, 3, 9, 13]])
exchangedMatrix = np.array([[1, 1, 1, 3],
                            [0, 1, 3, 4],
                            [0, 0, 1, 1]])
calculate = toolbox.gauss_elimination(matrix)
print("Result:")
print('Calculation result:\n', calculate)
print('Standard result\n', exchangedMatrix)


# 对高斯若尔当消元法的测试
print("\n")
print("###############################")
print("Test of Gauss Jordan Elimination algorithm.")
print("Result of P44 3nd:\n")
m1 = np.array([[-3, 8, 5, 1, 0, 0],
               [2, -7, 4, 0, 1, 0],
               [1, 9, -6, 0, 0, 1]])
print(toolbox.gauss_jordan_elimination(m1))

print("---------------------------")

print("Result of P44 7th:\n")
m2 = np.array([[2, -3, 1, 3],
               [-3, 5, -1, -6],
               [1, -1, 2, -2]])
result2 = toolbox.gauss_jordan_elimination(m2)
print("\tResult of 1st of P44 7th:\n")
print("\tMatrix of calculation")
for i in range(len(result2)):
    print("\t", result2[i])
print("\n")
print("\tIndependent variable:")
for i in range(len(result2)):
    print("\t", f"x{i + 1}: {result2[i][-1]}")

m3 = np.array([[2, -3, 1, 1, 0, 0],
               [-3, 5, -1, 0, 1, 0],
               [1, -1, 2, 0, 0, 1]])
result3 = toolbox.gauss_jordan_elimination(m3)
print("\n")
print("\tResult of 3nd of P44 7th:")
for i in range(len(result3)):
    print("\t", result3[i])


# 对追赶法的测试
print("\n")
print("###############################")
print("Test of chase algorithm.\n")
A = [-1, -1, -1, -1]
B = [4, 4, 4, 4, 4]
C = [-1, -1, -1, -1]
F = [100, 200, 200, 200, 100]
print("Result of P56 3rd:\n")
result = toolbox.chase(np.array(A), np.array(B), np.array(C), np.array(F))
print("Matrix of calculation by chase algorithm:")
print(result)
print("\n")
print("Independent variable:")
for i in range(len(result)):
    print(f"x{i + 1}: {result[i]}")


# 对道立特分解法的测试
m = np.array([[4, 1, -1, 0, 7],
              [1, 3, -1, 0, 8],
              [-1, -1, 5, 2, -4],
              [0, 0, 2, 4, 6]])
print("###############################")
print("Test of doolittle algorithm:\n")
print(toolbox.doolittle(m))


# 对高斯赛德尔迭代法的测试
A = np.array([[4, -1, 0, -1, 0, 0],
              [-1, 4, -1, 0, -1, 0],
              [0, -1, 4, 0, 0, -1],
              [-1, 0, 0, 4, -1, 0],
              [0, -1, 0, -1, 4, -1],
              [0, 0, -1, 0, -1, 4]])
f = np.array([0, 5, 0, 6, -2, 6])
e1 = 0.01
e2 = 0.00005
t = 1000
print("Test of G-S matrix form:")
print(toolbox.gauss_seidel_matrixForm(A, f, e2, t))
print("\n")

# 对松弛法的测试
A = np.array([[4, -1, 0],
              [-1, 4, -1],
              [0, -1, 4]])
b = np.array([1, 4, -3])
omega1 = 1.03
omega2 = 1
omega3 = 1.1
e = 0.000005
t = 1000
print("Test of SOR:")
print(toolbox.SOR(A, b, omega1, e, t))


# 对BM算法的测试
testBIn = "1 - x - math.sin(x)"
print("Test of BM algorithm:")
bm_result = toolbox.BM(testBIn, 0, 1, 0.001, 0.005, 10000)
print(f"最终结果：{bm_result[0]}")
print(f"迭代次数：{bm_result[1]}")


# 对试位法的测试
print("Test of false position algorithm:")
false_position_re = toolbox.false_position(testBIn, 0, 1, 0.01, 0.005, 10000)
print(f"最终结果：{false_position_re[0]}")
print(f"迭代次数：{false_position_re[1]}")


# 对BMZ算法的测试
testBIn1 = "6 * x ** 4 - 40 * x ** 2 + 9"
print("Test of BMZ algorithm:")
bmz_result = toolbox.BMZ(testBIn1, -5, 0.001, 0.1, 4, 0.1, 100)
print(f"所有根：{bmz_result[0]}")
print(f"迭代次数：{bmz_result[1]}")


#对简单迭代法的测试
testIter = "math.e ** (-x)"
print("Test of simple iteration algorithm:")
simpleIteration_result = toolbox.simple_iteration(testIter, 0.5, 0.0001, 0.0001, 10000)
print(f"最终结果：{simpleIteration_result[0]}")
print(f"迭代次数：{simpleIteration_result[1]}")


# 对埃特金法的测试
print("Test of steffensen algorithm:")
steffensen_result = toolbox.steffensen(testIter, 0.5, 0.0001, 0.0001, 10000)
print(f"最终结果：{steffensen_result[0]}")
print(f"迭代次数：{steffensen_result[1]}")


# 对牛顿法的测试
testNewton = "(x ** 2) ** (1 / 3)"
print("Test of Newton iteration algorithm:")
newtonSolution = toolbox.newton_iteration(testNewton, 1.5, 0.00001, 0.00001, 10000)
print(f"最终结果：{newtonSolution[0]}")
print(f"迭代次数：{newtonSolution[1]}")


# 对重根加速法的测试
print("Test of multiple solution algorithm:")
multi = toolbox.multiple_solutionNewton(testNewton, 1, 0.00001, 0.00001, 10000)
print(f"最终结果：{multi[0]}")
print(f"迭代次数：{multi[1]}")


# 对弦截法的测试
print("Test of secant algorithm:")
secantSolution = toolbox.secant(testNewton, 1, 0.00001, 0.00001, 10000)
print(f"最终结果：{secantSolution[0]}")
print(f"迭代次数：{secantSolution[1]}")


# 对牛顿下山法的测试
print("Test of Newton descent algorithm:")
des = toolbox.newton_descent(testNewton, 1.5, 0.00001, 0.00001, 0.001, 0.001, 10000)
print(f"最终结果：{des[0]}")
print(f"迭代次数：{des[1]}")


# 对插值法的测试
x_for_lagrange = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
y_for_lagrange = [0.4794, 0.6442, 0.7833, 0.8912, 0.9636, 0.9975, 0.9917, 0.9463]
lagrange = toolbox.Lagrange(x_for_lagrange, y_for_lagrange, 2)
print("拉格朗日插值法插值0.6：")
print(lagrange.interpolation(0.6))
print("拉格朗日插值法插值0.8：")
print(lagrange.interpolation(0.8))
print("拉格朗日插值法插值1：")
print(lagrange.interpolation(1))


# 对牛顿插值法的测试
x_for_newton = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
y_for_newton = [-0.916291, -0.693147, -0.510826, -0.357765, -0.223144, -0.105361]
newton = toolbox.Newton(x_for_newton, y_for_newton, 4)
node = 0.78
print("差分表：")
print("牛顿后插公式插值0.78：")
print(newton.interpolation(node))
print(newton._table)