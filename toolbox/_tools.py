import numpy as np


def my_max(matrix: np.ndarray) -> list and int:
    """
    输入矩阵，将
    """
    maxRow = 0
    maxColumn = 0
    maxElement = matrix[maxRow][maxColumn]
    for row in range(len(matrix)):
        for column in range(len(matrix[row]) - 1):
            if matrix[row][column] > maxElement:
                maxRow = row
                maxColumn = column
                maxElement = matrix[row][column]
    return [maxRow, maxColumn], maxElement


def float_transform(array_: np.ndarray, is_one_dimensional_array: bool) -> np.ndarray:
    """
    In order to compatible int and float, I use "float" type to change the type of all element in "array_" from
    int to float.
    """
    if is_one_dimensional_array:
        arr = [float(i) for i in array_]
    else:
        columns, rows = np.array(array_).shape
        floatArray = [float(i) for j in array_ for i in j]
        arr = np.array(floatArray).reshape(columns, rows)
    return arr
