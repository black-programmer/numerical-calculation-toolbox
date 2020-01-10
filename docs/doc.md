# 数值计算文档
## gauss_elimination
使用高斯消去法分解矩阵

例：
```python
import toolbox
import numpy as np
matrix = np.array([[1, 1, 1, 3],
                   [1, 2, 4, 7],
                   [1, 3, 9, 13]])
exchangedMatrix = np.array([[1, 1, 1, 3],
                            [0, 1, 3, 4],
                            [0, 0, 1, 1]])
calculate = toolbox.gauss_elimination(matrix)
print(calculate)
```
## doolittle
使用道立特分解法来分解矩阵

例：
```python
import toolbox
import numpy as np
m = np.array([[4, 1, -1, 0, 7],
              [1, 3, -1, 0, 8],
              [-1, -1, 5, 2, -4],
              [0, 0, 2, 4, 6]])
print(toolbox.doolittle(m))
```

## gauss_jordan_elimination
使用高斯若尔当消元法来分解矩阵

例：
```python
import toolbox
import numpy as np
m = np.array([[4, 1, -1, 0, 7],
              [1, 3, -1, 0, 8],
              [-1, -1, 5, 2, -4],
              [0, 0, 2, 4, 6]])
print("###############################")
print("Test of doolittle algorithm:\n")
print(toolbox.doolittle(m))
```

## chase
使用追赶法来分解矩阵

例：
```python
import toolbox
import numpy as np
A = [-1, -1, -1, -1]
B = [4, 4, 4, 4, 4]
C = [-1, -1, -1, -1]
F = [100, 200, 200, 200, 100]
result = toolbox.chase(np.array(A), np.array(B), np.array(C), np.array(F))
print(result)
```

## gauss_seidel_matrxForm
使用高斯赛德尔迭代法来解线性方程组

例：
```python
import toolbox
import numpy as np
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
print(toolbox.gauss_seidel_matrixForm(A, f, e2, t))
```

## SOR
使用松弛法来解线性方程组

例：
```python
import toolbox
import nu mpy
A = np.array([[4, -1, 0],
              [-1, 4, -1],
              [0, -1, 4]])
b = np.array([1, 4, -3])
omega1 = 1.03
omega2 = 1
omega3 = 1.1
e = 0.000005
t = 1000
print(toolbox.SOR(A, b, omega1, e, t))
```

## BM
使用BM算法来解非线性方程

例：
```python
import toolbox
import numpy as np
testBIn = "1 - x - math.sin(x)"
print("Test of BM algorithm:")
bm_result = toolbox.BM(testBIn, 0, 1, 0.001, 0.005, 10000)
print(bm_result)
```

## BMZ
使用BMZ算法来解非线性方程

例：
```python
import toolbox
import numpy as np
testBIn1 = "6 * x ** 4 - 40 * x ** 2 + 9"
print("Test of BMZ algorithm:")
bmz_result = toolbox.BMZ(testBIn1, -5, 0.001, 0.1, 4, 0.1, 100)
print(bmz_result)
```

## simple_iteration
使用简单迭代法解非线性方程组

例：
```python
import toolbox
import numpy as np
testIter = "math.e ** (-x)"
print("Test of simple iteration algorithm:")
simpleIteration_result = toolbox.simple_iteration(testIter, 0.5, 0.0001, 0.0001, 10000)
print(simpleIteration_result)
```

## steffense
使用埃特金法解非线性方程组

```python
import toolbox
import numpy as np
steffensen_result = toolbox.steffensen(testIter, 0.5, 0.0001, 0.0001, 10000)
print(steffensen_result)
```

## newton_iteration
使用牛顿法解非线性方程组

例：
```python
import toolbox
import numpy as np
testNewton = "(x ** 2) ** (1 / 3)"
newtonSolution = toolbox.newton_iteration(testNewton, 1.5, 0.00001, 0.00001, 10000)
print(nuewtonSolution)
```

## multiple_solutionNewton
使用重根加速法解非线性方程组

例：
```python
import toolbox
import numpy as np
multi = toolbox.multiple_solutionNewton(testNewton, 1, 0.00001, 0.00001, 10000)
```

## secant
使用弦截法解非线性方程组

例：
```python
import toolbox
import numpy as np
secantSolution = toolbox.secant(testNewton, 1, 0.00001, 0.00001, 10000)
print(secantSolution)
```

## newton_descent
使用牛顿下山法解非线性方程组

例：
```python
import toolbox
import numpy as np
des = toolbox.newton_descent(testNewton, 1.5, 0.00001, 0.00001, 0.001, 0.001, 10000)
print(des)
```

## Lagrange
拉格朗日插值法

例：
```python
import toolbox
import numpy as np
x_for_lagrange = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
y_for_lagrange = [0.4794, 0.6442, 0.7833, 0.8912, 0.9636, 0.9975, 0.9917, 0.9463]
lagrange = toolbox.Lagrange(x_for_lagrange, y_for_lagrange, 2)
print(lagrange.interpolation(0.6))
print(lagrange.interpolation(0.8))
print(lagrange.interpolation(1))
```

## Newton
牛顿插值法

例：
```python
import toolbox
import numpy as np
x_for_newton = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
y_for_newton = [-0.916291, -0.693147, -0.510826, -0.357765, -0.223144, -0.105361]
newton = toolbox.Newton(x_for_newton, y_for_newton, 4)
node = 0.78
print(newton.interpolation(node))
```