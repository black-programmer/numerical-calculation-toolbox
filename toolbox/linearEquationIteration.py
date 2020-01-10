from ._tools import *


def gauss_seidel_matrixForm(arr: np.ndarray, b: np.ndarray, epsilon: float, time: int) -> tuple:
    """
    This method can solute a linear equation by Gauss-Seider algorithm.

    Example:
            A = np.array([[4, -1, 0, -1, 0, 0],
                          [-1, 4, -1, 0, -1, 0],
                          [0, -1, 4, 0, 0, -1],
                          [-1, 0, 0, 4, -1, 0],
                          [0, -1, 0, -1, 4, -1],
                          [0, 0, -1, 0, -1, 4]])
            f = np.array([0, 5, 0, 6, -2, 6])
            e = 0.00005
            t = 10
            result = gaussSeidel(A, f, e, t)
            print(result)

    :param arr: np.ndarray type, the coefficient matrix of objective linear equation.
    :param b: np.ndarray type, the constant sequence of objective linear equation
    :param epsilon: float type, the accuracy.
    :param time: int type, the max iteration time.
    :return: tuple type, the first element is the result of objective linear equation
        and the second element is the iteration time of this calculation.
    """

    arr_ = float_transform(arr, False)
    # Using "floatTransform" in "tools" module to transform the each element type of
    # "arr" to float type.

    x = [0] * len(arr)
    # The initial value of each independent variable of objective equation.

    sigmaA = []
    k = 1

    while k <= time:
        xi = np.array([ele for ele in x])
        # If using "xi = x" here, the value of "err" will be zero forever, because
        # the assignment operation in Python is "appoint transfer" not "value transfer",
        # That means "xi" and "x" are both point the same value.

        for i in range(len(arr_)):
            for j in range(len(arr_[0])):
                if i != j:
                    sigmaA.append(arr_[i][j] * x[j])
            sumSigmaA = sum(sigmaA)
            sigmaA = []
            x[i] = (b[i] - sumSigmaA) / arr_[i][i]

        err = max(abs(np.array(xi) - np.array(x)))
        print(x)
        if err < epsilon:
            break

        k += 1

    return x, k


def SOR(arr: np.ndarray, b: np.ndarray, omega: float, epsilon: float, time: int) -> tuple:
    x = [0] * len(arr)
    # The initial value of each independent variable of objective equation.

    sigmaA = []
    k = 1

    while k <= time:
        xi = [ele for ele in x]
        # If using "xi = x" here, the value of "err" will be zero forever, because
        # the assignment operation in Python is "appoint transfer" not "value transfer",
        # That means "xi" and "x" are both point the same value.

        for i in range(len(arr)):
            for j in range(len(arr[0])):
                if i != j:
                    sigmaA.append(arr[i][j] * x[j])
            sumSigmaA = sum(sigmaA)
            sigmaA = []
            x[i] = omega * (b[i] - sumSigmaA) / arr[i][i]

        err = max(abs(np.array(xi) - np.array(x)))
        print(x)
        if err < epsilon:
            break

        k += 1

    return x, k


if __name__ == '__main__':
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
    print(gauss_seidel_matrixForm(A, f, e2, t))
    print("\n")

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
    print(SOR(A, b, omega1, e, t))
