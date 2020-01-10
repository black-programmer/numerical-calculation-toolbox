import math
from functools import partial
import numpy as np
import matplotlib.pyplot as plt
from scipy.misc import derivative


def _function(fun: str, x: int or float):
    """
    Transform the string to executable form.

    :param fun: a one parameter function, it must be a string.
    :param x: the parameter of the function..
    :return: the executable form of the string function.
    """
    return eval(fun)


def BM(
        fun: str,
        a: int or float,
        b: int or float,
        eps1: float,
        eps2: float,
        max_time: int
) -> tuple:
    """
    BM algorithm.

    :param fun: the function which will be solved by BM, it must be a string.
    :param a: the left boundary you assign.
    :param b: the right boundary you assign.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is a list which
        includes the single solution of objective equation and the second
        element is the iteration time.
    """
    f = partial(_function, fun)
    # The partial function, it can fix the parameters for a objective function.

    time = 1
    while time < max_time:
        x = (a + b) / 2
        xe = (b - a) / 2
        print(f"第{time}次迭代:{x}")

        if abs(xe) < eps2 and abs(f(x)) < eps1:
            return x, time
        else:
            if f(a) * f(x) < 0:
                b = x
            else:
                a = x
        time += 1


def BMZ(
        fun: str,
        a: int or float,
        eps1: float,
        eps2: float,
        k: int,
        dx: int or float,
        max_time: int
) -> tuple:
    """
    BMZ algorithm.

    :param fun: the function which will be solved by BM, it must be a string.
    :param a: the left boundary you assign.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param k: the number of all solutions.
    :param dx: the step length.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is a list which
        includes the single solution of objective equation and the second
        element is the iteration time.
    """
    f = partial(_function, fun)
    # The partial function, it can fix the parameters for a objective function.

    all_solution = []

    time = 1
    while time < max_time:
        x = a + time * dx

        if abs(f(x) - 0) < eps1:
            all_solution.append(x)

        elif f(x) * f(x + dx) > 0:
            # The judgement condition in book is "f(x) * f(x + 1) > 0", but
            # while I debugged here, I discovered it will ignore one solution,
            # so I had to modify it.

            time += 1
            pass
        else:
            try:
                all_solution.append(
                    BM(fun, x, a + (time + 1) * dx, eps1, eps2, max_time)[0]
                )
                time += 1
            except Exception:
                time += 1
                pass

        if len(all_solution) == k:
            break

    return all_solution, time


def false_position(
        fun: str,
        a: int or float,
        b: int or float,
        eps1: float,
        eps2: float,
        max_time: int
) -> tuple:
    """
    Using false position algorithm to solve the function.

    :param fun: the function which will be solved by false position
        algorithm, it must be a string.
    :param a: the left boundary.
    :param b: the right boundary.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: the single solution of the function.
    """

    def _getX1(fun, x_left, x_right):
        x = x_left - (
                (x_right - x_left) / (fun(x_right) - fun(x_left))
            ) * fun(x_left)
        return x

    f = partial(_function, fun)
    # The partial function, it can fix the parameters for a objective function.

    getX1 = partial(_getX1, f)
    x1 = getX1(a, b)
    x = a
    xe = b

    time = 1
    while time < max_time:
        if x * x1 < 0:
            xe = x1
            x1 = getX1(x, xe)
        else:
            x = x1
            x1 = getX1(x, xe)

        if abs((xe - x) / 2) < eps1 or abs(f(x)) < eps2:
            break
        time += 1

    return x, time


def simple_iteration(
        fun: str,
        initial_value: int or float,
        eps1: int or float,
        eps2: int or float,
        max_time: int
) -> tuple:
    """
    Using simple iteration algorithm to calculate the solution of nonlinear
    equation which you assign. There is a little inconvenience, the parameter
    "fun" must be the iteration form but not the objective function. Be example,
    if you want to get the solution of "x ** 2", you must make the "fun" is
    iteration form like "x ** 2 + x" or other iteration forms.

    :param fun: The iteration form of your objective function, it must be str
        type.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is a list which include
        all solutions of the objective equation and the second element is the
        iteration time.
    """

    f = partial(_function, fun)
    x = initial_value
    solution_sequence = [x]

    time = 1
    while time < max_time:
        x = f(x)
        solution_sequence.append(x)

        if abs(solution_sequence[-1] - solution_sequence[-2]) < eps1 \
                or abs(f(solution_sequence[-1])) < eps2:
            break

        time += 1
    return solution_sequence, time


def steffensen(
        fun: str,
        initial_value: float or int,
        eps1: float or int,
        eps2: float or int,
        max_time: int
) -> tuple:
    """
    This function implements Steffensen algorithm which can solve a nonlinear
    equation, compare to simple iteration algorithm, it has more fast convergence
    rate. And, same as simple iteration, you must transform the objective equation
    to iteration form.

    :param fun: the iteration form of the objective equation.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is the solution and the
        second element is the iteration time.
    """
    f = partial(_function, fun)
    x = [initial_value]
    print(f"第1次：{x[-1]}， f(x) = {f(x[-1]) - x[-1]}")

    time = 1
    while time < max_time:
        x1 = f(x[-1])
        x2 = f(x1)
        x.append(x[-1] - ((x1 - x[-1]) ** 2) / (x2 - 2 * x1 + x[-1]))

        if abs(x[-1] - x[-2]) < eps1 or abs(f(x[-1])) < eps2:
            break

        print(f"第{time}次：{x[-1]}， f(x) = {f(x[-1]) - x[-1]}")

        time += 1
    print(f"最后一次：{x[-1]}，f(x) = {f(x[-1]) - x[-1]}")

    return x, time


def newton_iteration(
        fun: str,
        initial_value: int or float,
        eps1: int or float,
        eps2: int or float,
        max_time: int
) -> tuple:
    """
    This function is a implementation of Newton iteration algorithm. The
    parameter "fun" is the iteration form of the objective equation.

    :param fun: the iteration form of the objective equation.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is the solution and the
        second element is the iteration time.
    """
    f = partial(_function, fun)
    x = initial_value
    print(f"第1次迭代结果：{x}，f(x) = {f(x)}")

    time = 1
    while time < max_time:
        xOld = x

        fd = (f(x + 0.0000001) - f(x)) / 0.0000001
        # Using the definition to calculate the value of derivatives of the
        # objective equation.

        x -= f(x) / fd

        if abs(x - xOld) < eps1 and abs(f(x)) < eps2:
            break
        print(f"第{time}次迭代结果：{x}，f(x) = {f(x)}")

        time += 1

    return x, time


def multiple_solutionNewton(
        fun: str,
        initial_value: int or float,
        eps1: int or float,
        eps2: int or float,
        max_time: int
) -> tuple:
    """
    The algorithm of this function implementation is multiple root acceleration.

    :param fun: the iteration form of the objective equation.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is the solution and the
        second element is the iteration time.
    """
    f = partial(_function, fun)
    x = initial_value

    time = 1
    while time < max_time:
        xOld = x

        fd = derivative(f, x, dx=0.0001, n=1)
        fdd = derivative(f, x, dx=0.00001, n=2)
        # In order to calculate the value of two derivative of the objective
        # equation, I have to use the function "derivative" in "scipy". The
        # variable "fd" is the value of one derivative of the objective and
        # the "fdd" is the value of two derivative.

        x -= (f(x) * fd) / (fd ** 2 - f(x) * fdd)

        if abs(x - xOld) < eps1 and abs(f(x) < eps2):
            break

        time += 1

    return x, time


def secant(
        fun: str,
        initial_value: int or float,
        eps1: int or float,
        eps2: int or float,
        max_time: int
) -> tuple:
    """
    The algorithm this function implementation is the improvement of Newton
    iteration algorithm, this algorithm use difference quotient to replace
    differential quotient.


    :param fun: the iteration form of the objective equation.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is the solution and the
        second element is the iteration time.
    """
    f = partial(_function, fun)
    x = [initial_value]
    x.append(x[0] - f(x[0]) / (((f(x[0]) + 0.0000001) - f(x[0])) / 0.0000001))

    time = 1
    while time < max_time:
        xNew = x[-1] - (f(x[-1]) / (f(x[-1]) - f(x[-2]))) * (x[-1] - x[-2])
        x.append(xNew)

        if abs(x[-1] - x[-2]) < eps1 and abs(f(x[-1])) < eps2:
            break

        time += 1

    return x[-1], time


def newton_descent(
        fun: str,
        initial_value: int or float,
        eps1: int or float,
        eps2: int or float,
        eps3: int or float,
        eps4: int or float,
        max_time: int
) -> tuple:
    """
    This function implements the improvement of Newton iteration, it has more fast
    convergence speed.

    :param fun: the iteration form of the objective equation.
    :param initial_value: the initial value.
    :param eps1: the minimum boundary of the solution.
    :param eps2: the minimum boundary of the function's value of solution.
    :param eps3: the minimum value of t.
    :param eps4: the correction of initial value.
    :param max_time: the max iteration time.
    :return: a tuple with two elements, the first element is the solution and the
        second element is the iteration time.
    """
    f = partial(_function, fun)
    x = [initial_value]
    t = 1

    time = 1
    while time < max_time:
        fd = (f(x[-1] + 0.0000001) - f(x[-1])) / 0.0000001
        x_ = x[-1] - t * (f(x[-1]) / fd)

        if abs(f(x_) - f(x[-1])) < eps1 or abs(f(x_)) < eps2:
            return x_, time

        if abs(f(x_)) < abs(f(x[-1])):
            x.append(x_)
            print(f"第{time}次迭代结果：{x[-1]}")
            time += 1
        else:
            if t > eps3:
                t /= 2
                time += 1
                continue
            else:
                x.append(x_ + eps4)
                print(f"第{time}次迭代结果：{x[-1]}")
                t /= 2
                time += 1

    return x[-1], time


if __name__ == '__main__':
    testBIn = "1 - x - math.sin(x)"
    print("Test of BM algorithm:")
    bm_result = BM(testBIn, 0, 1, 0.001, 0.005, 10000)
    print(f"最终结果：{bm_result[0]}")
    print(f"迭代次数：{bm_result[1]}")

    print("Test of false position algorithm:")
    false_position_re = false_position(testBIn, 0, 1, 0.01, 0.005, 10000)
    print(f"最终结果：{false_position_re[0]}")
    print(f"迭代次数：{false_position_re[1]}")

    testBIn1 = "6 * x ** 4 - 40 * x ** 2 + 9"
    print("Test of BMZ algorithm:")
    bmz_result = BMZ(testBIn1, -5, 0.001, 0.1, 4, 0.1, 100)
    print(f"所有根：{bmz_result[0]}")
    print(f"迭代次数：{bmz_result[1]}")

    testIter = "math.e ** (-x)"
    print("Test of simple iteration algorithm:")
    simpleIteration_result = simple_iteration(testIter, 0.5, 0.0001, 0.0001, 10000)
    print(f"最终结果：{simpleIteration_result[0]}")
    print(f"迭代次数：{simpleIteration_result[1]}")
    print("Test of steffensen algorithm:")
    steffensen_result = steffensen(testIter, 0.5, 0.0001, 0.0001, 10000)
    print(f"最终结果：{steffensen_result[0]}")
    print(f"迭代次数：{steffensen_result[1]}")
    plt.plot(simpleIteration_result[0], label="simple iteration")
    plt.plot(steffensen_result[0], label="steffensen")
    plt.show()

    # testNewton = "x - (1 / 2) * (x + 3 / x)"  # 3 ^ (1 / 2) iteration form.
    testNewton = "(x ** 2) ** (1 / 3)"
    # test = "(x - 1.5) ** 2"
    print("Test of Newton iteration algorithm:")
    newtonSolution = newton_iteration(testNewton, 1.5, 0.00001, 0.00001, 10000)
    print(f"最终结果：{newtonSolution[0]}")
    print(f"迭代次数：{newtonSolution[1]}")
    print("Test of multiple solution algorithm:")
    multi = multiple_solutionNewton(testNewton, 1, 0.00001, 0.00001, 10000)
    print(f"最终结果：{multi[0]}")
    print(f"迭代次数：{multi[1]}")

    print("Test of secant algorithm:")
    secantSolution = secant(testNewton, 1, 0.00001, 0.00001, 10000)
    print(f"最终结果：{secantSolution[0]}")
    print(f"迭代次数：{secantSolution[1]}")
    print("Test of Newton descent algorithm:")
    des = newton_descent(testNewton, 1.5, 0.00001, 0.00001, 0.001, 0.001, 10000)
    print(f"最终结果：{des[0]}")
    print(f"迭代次数：{des[1]}")
    #
    # t = "(x ** 2) ** (1 / 3)"
    # x = np.linspace(-10, 10, 10000)
    # plt.plot(x, list(map(partial(_function, t), x)))
    # plt.show()

    # def e(x):
    #     return math.e ** -x
    # def f(x):
    #     ff = x - (e(x) - x) ** 2 / (e(e(x)) - 2 * e(x) + x)
    #     return ff
    # print(derivative(f, 0.5671190400572149, dx=0.0001, n=1))
    # print(derivative(f, 0.5671190400572149, dx=0.0001, n=2))
    # def ff(x):
    #     return abs(derivative(e, x, dx=0.0001, n=1))
    #
    # xx = np.linspace(0, 1, 1000)
    #
    # yy = max(list(map(ff, xx)))
    # print(yy)

    # x = np.linspace(0, 1, 10000)
    # plt.plot(x, list(map(f, x)))
    # plt.show()
    # print(e(0), e(1))
    # print(f(1))
