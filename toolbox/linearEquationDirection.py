from ._tools import *


class ElementaryTransformation:
    """
    In order to use conveniently, I putted three basic elementary transformations in this class and make them
    as static method, that represent I can use them without instantiation.
    """

    @staticmethod
    def multiple(num: int, array_: np.ndarray) -> np.ndarray:
        arr = np.array(array_)
        multipleArray = num * arr
        return multipleArray

    @staticmethod
    def multiple_add(num: int, multiple_index: int, add_index: int, array_: np.ndarray) -> np.ndarray:
        if multiple_index != add_index:
            arr = np.array(array_)
            newLine = (num * arr[multiple_index]) + array_[add_index]
            arr = np.delete(arr, add_index, 0)
            arr = np.insert(arr, add_index, newLine, 0)
            return arr
        else:
            return array_

    @staticmethod
    def exchange_rows(change_index: int, changed_index: int, array_: np.ndarray) -> np.ndarray:
        arr = list(array_)
        arr[changed_index], arr[change_index] = arr[change_index], arr[changed_index]
        return np.array(arr)


def gauss_elimination(array_: np.ndarray) -> np.ndarray:
    arr = float_transform(array_, False)

    def elimination(arr_=arr, i=0):
        for k in range(i, len(arr_) - 1):
            if arr_[k][k] == 0 and arr_[k + 1][k] != 0:

                exchangeArr = ElementaryTransformation.exchange_rows(k, len(arr_) - 1, arr_)
                elimination(exchangeArr, k)
                # if a[k][k] = 0 and [k+1][k] != 0, then these to rows will be exchanged. The argument will insure this
                # function continue to run.

            elif arr_[k][k] == 0 and arr_[k + 1][k] == 0:
                if k + 1 < len(arr_) - 1:
                    elimination(arr_, k + 1)
                else:
                    return arr_

            for i in range(k + 1, len(arr_)):
                arr_[i][k] = arr_[i][k] / arr_[k][k]
                for j in range(k + 1, len(arr_[k])):
                    arr_[i][j] -= arr_[i][k] * arr_[k][j]
                    print('Each result of calculation:\n', arr_, '\n')
        return arr_

    return elimination()


def gauss_jordan_elimination(arr_: np.ndarray) -> np.ndarray:
    """
    This function is based on Gauss-Jordan elimination algorithm:
        k = 0, 1, 2, ..., n - 1
        > a[k][j] <- a[k][j] / a[k][k] (j = k, k + 1, ..., n)
        > i != k
        > a[i][j] <- a[i][j] - a[i][k] * a[k][j]

    :param arr_: A np.ndarray type matrix, which is the augmented matrix of objective matrix.
    :return: A np.array type array, it have two parts and the first part is result matrix,
        the second part is index of independent variable "x" of objective function.
    """
    arr = float_transform(arr_, False)
    for k in range(len(arr)):
        for j in range(k, len(arr[0])):
            arr[k][j: len(arr[0])] /= arr[k][k]
            for i in range(len(arr)):
                if i != k:
                    arr[i][j: len(arr[0])] -= arr[i][k] * arr[k][j: len(arr[0])]

    return arr


def doolittle(arr: np.ndarray, f: np.ndarray = None, choose: bool = False) -> np.ndarray:
    """
    This function based on Doolittle algorithm can solute some equation or do L-U decomposition
    for the matrix you defined.

    :param arr: np.ndarray type, the coefficient matrix for the objective equation.
    :param f: np.ndarray type, the constant sequence for the objective equation.
    :param choose: bool type, if you want to get the result of the objective equation,
        you must make it be "True", or you can ignore it.
    :return: np.ndarray type, the result of L-U decomposition or solution of the equation.
    """
    arr_ = float_transform(arr, False)
    if f is not None:
        f_ = float_transform(f, True)

    columns = []
    rows = []
    for r in range(len(arr_)):
        for j in range(r, len(arr_[0])):
            for k in range(r):
                columns.append(arr_[r][k] * arr_[k][j])
            sumC = sum(columns)
            columns = []
            arr_[r][j] -= sumC
        for i in range((r + 1), len(arr_)):
            for k in range(r):
                rows.append(arr_[i][k] * arr_[k][r])
            sumR = sum(rows)
            rows = []
            arr_[i][r] = (arr_[i][r] - sumR) / arr_[r][r]

    if f is not None:
        AY = []
        for i in range(len(f_)):
            for k in range(i):
                AY.append(arr_[i][k] * f_[k])
            sumAY = sum(AY)
            AY = []
            f_[i] -= sumAY

        AX = []
        for i in range((len(f_) - 1), -1, -1):
            for k in range((i + 1), len(arr_)):
                AX.append(arr_[i][k] * f_[k])
            sumAX = sum(AX)
            AX = []
            f_[i] = (f_[i] - sumAX) / arr_[i][i]

    if choose is True:
        return f_
    else:
        return arr_


def chase(a: np.ndarray, b: np.ndarray, c: np.ndarray, f: np.ndarray) -> list:
    """
    This method can solve tridiagonal equation by pursuit algorithm:
        first, I must use Crout decompose algorithm to get "L" and "U":
            gamma[i] = a[i] (i = 1, 2, ..., n - 1)
            alpha[0] = b[0], beta[1] = c[0] / alpha[0]
            > alpha[i] = b[i] - gamma[i] * beta[i - 1]
            > beta[i] = c[i] / alpha[i] (i = 1, 2, ..., n - 2)
            alpha[-1] = b[-1] - gamma[-1] * beta[-2]
        Then use pursuit algorithm:
            y[1] = f[1] / alpha[1], y[i] = (f[i] - gamma[i - 1]) / alpha[i] (i = 1, 2, ..., n)
            x[-1] = y[-1], x[i] = y[i] - beta[i] * x[i + 1] (i = n - 2, n - 3, ..., 2, 1, 0)

    :param a: this parameter correspond to left element of diagonal element of objective function,
        it should be a np.ndarray type array.
    :param b: this parameter correspond to diagonal element of objective function, it should be a
        np.ndarray type array.
    :param c: this parameter correspond to diagonal element of objective function, it should be a
        np.ndarray type array.
    :param f: this parameter correspond to constant sequence of objective function, it should be
        a np.ndarray type array.
    :return: the result of objective function, it is list type.
    """

    gamma = a
    alpha = [b[0]]
    beta = [c[0] / alpha[0]]

    for i in range(1, len(b) - 1):
        # Append element from second location, so starting from 1.
        alpha.append(b[i] - gamma[i - 1] * beta[i - 1])
        beta.append(c[i] / alpha[i])
        # gamma[i] = a[i] (i = 1, 2, ..., n - 1)
        # alpha[0] = b[0], beta[1] = c[0] / alpha[0]
    alpha.append(b[-1] - gamma[-1] * beta[-2])
    # alpha[-1] = b[-1] - gamma[-1] * beta[-2], the index "-1" represent the last element of list alpha.

    y = [f[0] / alpha[0]]
    for i in range(1, len(f)):
        y.append((f[i] - (gamma[i - 2] * y[i - 1])) / alpha[i])
    x = [0, 0, 0, 0, y[-1]]
    # Initialize the list "x", because this step appends the element is back to front, so we must modify
    # list "x" from the last element to first element and not use "append".
    for i in range(len(f) - 2, -1, -1):
        print(x)
        x[i] = y[i] - beta[i] * x[i + 1]

    return x


if __name__ == '__main__':
    print("\n")
    print("Test of Gauss Elimination algorithm.")
    matrix = np.array([[1, 1, 1, 3],
                       [1, 2, 4, 7],
                       [1, 3, 9, 13]])
    exchangedMatrix = np.array([[1, 1, 1, 3],
                                [0, 1, 3, 4],
                                [0, 0, 1, 1]])
    calculate = gauss_elimination(matrix)
    print("Result:")
    print('Calculation result:\n', calculate)
    print('Standard result\n', exchangedMatrix)

    print("\n")
    print("###############################")
    print("Test of Gauss Jordan Elimination algorithm.")
    print("Result of P44 3nd:\n")
    m1 = np.array([[-3, 8, 5, 1, 0, 0],
                   [2, -7, 4, 0, 1, 0],
                   [1, 9, -6, 0, 0, 1]])
    print(gauss_jordan_elimination(m1))

    print("---------------------------")

    print("Result of P44 7th:\n")
    m2 = np.array([[2, -3, 1, 3],
                   [-3, 5, -1, -6],
                   [1, -1, 2, -2]])
    result2 = gauss_jordan_elimination(m2)
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
    result3 = gauss_jordan_elimination(m3)
    print("\n")
    print("\tResult of 3nd of P44 7th:")
    for i in range(len(result3)):
        print("\t", result3[i])

    print("\n")
    print("###############################")
    print("Test of chase algorithm.\n")
    A = [-1, -1, -1, -1]
    B = [4, 4, 4, 4, 4]
    C = [-1, -1, -1, -1]
    F = [100, 200, 200, 200, 100]
    print("Result of P56 3rd:\n")
    result = chase(np.array(A), np.array(B), np.array(C), np.array(F))
    print("Matrix of calculation by chase algorithm:")
    print(result)
    print("\n")
    print("Independent variable:")
    for i in range(len(result)):
        print(f"x{i + 1}: {result[i]}")

    m = np.array([[4, 1, -1, 0, 7],
                  [1, 3, -1, 0, 8],
                  [-1, -1, 5, 2, -4],
                  [0, 0, 2, 4, 6]])
    print("###############################")
    print("Test of doolittle algorithm:\n")
    print(doolittle(m))
