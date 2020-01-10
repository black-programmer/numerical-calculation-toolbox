import math
from abc import ABC, abstractmethod


class InterpolationBase(ABC):
    def __init__(self, x_list: list, y_list: list, n_: int):
        self.x_list = x_list
        self.y_list = y_list
        if n_ == 0:
            raise ValueError(
                "the number of interpolation can't be zero."
            )
        if n_ > len(self.x_list) - 1:
            raise ValueError(
                "the number of interpolation must " +
                "smaller than length of x_list subtract one."
            )
        self.n = n_

    def _check_index(self, x: int or float) -> int:
        for i in range(len(self.x_list) - 1):
            if self.x_list[i] < x < self.x_list[i + 1]:
                index = i
                break
            elif self.x_list[i] == x:
                return i
            elif self.x_list[i + 1] == x:
                return i + 1
        else:
            raise ValueError(
                "Value of node must greater than the minimum " +
                "value and smaller than the maximum value of \"x_list\"."
            )

        return index

    @abstractmethod
    def interpolation(self, x: int or float) -> int or float:
        ...


class Lagrange(InterpolationBase):
    def __init__(self, x_list: list, y_list: list, n_: int):
        super().__init__(x_list, y_list, n_)
        self._x_node_list = []
        self._y_node_list = []

    def _create_interpolation_base_function(
            self, x_list: list, y_list: list
    ) -> str:
        numerator = []
        denominator = []
        function_part = []
        for each in x_list:
            numerator.append(f"(x - {each})")
            num = [ele for ele in x_list if ele != each]
            d = f"({each} - {num[0]})"
            for i in range(1, len(num)):
                d += f" * ({each} - {num[i]})"
            denominator.append(d)
        for i in range(len(x_list)):
            f1 = [ele for ele in numerator if ele != numerator[i]]
            func = f1[0]
            for index in range(1, len(f1)):
                func += f" * {f1[index]}"
            func += f" / ({denominator[i]})"
            func = f"{y_list[i]} * {func}"
            function_part.append(func)
        f = function_part[0]
        for i in range(1, len(function_part)):
            f += f" + {function_part[i]}"

        return f

    def _create_interpolation_function(
            self, x: int or float
    ) -> str or int or float:
        index = super()._check_index(x)
        if index < self.n:
            self._x_node_list = self.x_list[: self.n + 1]
            self._y_node_list = self.y_list[: self.n + 1]
        elif index + 1 > len(self.x_list) - self.n:
            self._x_node_list = self.x_list[-self.n - 1:]
            self._y_node_list = self.y_list[-self.n - 1:]
        else:
            start = math.ceil(math.sqrt(index))
            self._x_node_list = self.x_list[start: start + self.n + 1]
            self._y_node_list = self.y_list[start: start + self.n + 1]

        return self._create_interpolation_base_function(
            self._x_node_list, self._y_node_list
        )

    def _eval_function(self, func: str, x: float or int):
        return eval(func)

    def interpolation(self, x: int or float) -> int or float:
        result = self._create_interpolation_function(x)
        print(result)
        if isinstance(result, float or int):
            return result
        else:
            return self._eval_function(result, x)


class Newton(InterpolationBase):
    def __init__(self, x_list: list, y_list: list, n_: int):
        super().__init__(x_list, y_list, n_)
        self._table = self._create_difference_table()

    def _create_difference_table(self) -> dict:
        difference_table = {"fx": self.y_list}
        difference_list = []
        for i in range(len(difference_table["fx"]) - 1):
            difference_list.append(
                difference_table["fx"][i + 1] - difference_table["fx"][i]
            )
        difference_table["level1"] = difference_list
        for i in range(2, self.n + 1):
            difference_list = []
            for j in range(len(difference_table[f"level{i - 1}"]) - 1):
                difference_list.append(
                    difference_table[f"level{i - 1}"][j + 1] -
                    difference_table[f"level{i - 1}"][j]
                )
            difference_table[f"level{i}"] = difference_list

        return difference_table

    def _forward_interpolation(self, t: int or float):
        difference_coefficient = []
        for each in list(self._table.items()):
            difference_coefficient.append(each[1][0])
        value = difference_coefficient[0]
        for i in range(1, len(difference_coefficient)):
            for j in range(1, i):
                t *= t - j
            factorial_ = math.factorial(i)
            value += (t / factorial_) * difference_coefficient[i]
            t = 1

        return value

    def _backward_interpolation(self, t: int or float):
        difference_coefficient = []
        for each in list(self._table.items()):
            difference_coefficient.append(each[1][-1])
        value = difference_coefficient[0]
        for i in range(1, len(difference_coefficient)):
            for j in range(1, i):
                t *= t + j
            factorial_ = math.factorial(i)
            value += (t / factorial_) * difference_coefficient[i]
            t = 1

        return value

    def interpolation(self, x: int or float) -> int or float:
        index = super()._check_index(x)
        index_mid = len(self.x_list) >> 1
        if index < index_mid - 1:
            t = (x - self.x_list[0]) / (self.x_list[1] - self.x_list[0])
            value = self._forward_interpolation(t)
        else:
            t = (x - self.x_list[-1]) / (self.x_list[-1] - self.x_list[-2])
            value = self._backward_interpolation(t)

        return value


if __name__ == '__main__':
    x_for_lagrange = [0.5, 0.7, 0.9, 1.1, 1.3, 1.5, 1.7, 1.9]
    y_for_lagrange = [0.4794, 0.6442, 0.7833, 0.8912, 0.9636, 0.9975, 0.9917, 0.9463]
    lagrange = Lagrange(x_for_lagrange, y_for_lagrange, 2)
    print("拉格朗日插值法插值0.6：")
    print(lagrange.interpolation(0.6))
    print("拉格朗日插值法插值0.8：")
    print(lagrange.interpolation(0.8))
    print("拉格朗日插值法插值1：")
    print(lagrange.interpolation(1))

    x_for_newton = [0.4, 0.5, 0.6, 0.7, 0.8, 0.9]
    y_for_newton = [-0.916291, -0.693147, -0.510826, -0.357765, -0.223144, -0.105361]
    newton = Newton(x_for_newton, y_for_newton, 4)
    node = 0.78
    print("差分表：")
    print("牛顿后插公式插值0.78：")
    print(newton.interpolation(node))
    print(newton._table)