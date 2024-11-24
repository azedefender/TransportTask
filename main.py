import numpy as np

# Определяем матрицу затрат
c = np.array([
    [40, 30, 50, 45, 70, 0],
    [42, 28, 70, 44, 60, 0],
    [20, 36, 62, 71, 44, 0],
    [64, 41, 21, 37, 56, 0],
    [34, 38, 52, 61, 48, 0],
    [42, 27, 64, 48, 58, 0],
])

# Определяем векторы предложения и спроса
a = np.array([80, 160, 200, 70, 110, 1530])
b = np.array([80, 120, 200, 70, 110, 750])

# Создаем начальное решение с помощью метода северо-западного угла
def north_west_corner_rule(a, b, c):
    m, n = len(a), len(b)
    x = np.zeros((m, n))

    i, j = 0, 0
    while i < m and j < n:
        if a[i] < b[j]:
            x[i][j] = a[i]
            b[j] -= a[i]
            i += 1
        else:
            x[i][j] = b[j]
            a[i] -= b[j]
            j += 1

    return x

# Создаем начальное допустимое решение
initial_solution = north_west_corner_rule(a.copy(), b.copy(), c)

# Рассчитываем общую стоимость начального решения
total_cost = np.sum(initial_solution * c[:initial_solution.shape[0], :initial_solution.shape[1]])

# Функция для расчета потенциалов
def calculate_potentials(x, c):
    m, n = x.shape
    u = np.zeros(m)
    v = np.zeros(n)
    u[0] = 0  # Устанавливаем u1 = 0

    # Шаг 1: Рассчитываем v на основе u
    for i in range(m):
        for j in range(n):
            if x[i][j] > 0:
                v[j] = c[i][j] - u[i]

    # Шаг 2: Рассчитываем u на основе v
    for j in range(n):
        for i in range(m):
            if x[i][j] > 0:
                u[i] = c[i][j] - v[j]

    return u, v

# Функция для проверки оптимальности
def check_optimality(x, c):
    m, n = x.shape
    u, v = calculate_potentials(x, c)

    for i in range(m):
        for j in range(n):
            if x[i][j] == 0:  # Проверяем только не основные переменные
                if u[i] + v[j] > c[i][j]:
                    return False  # Не оптимально
    return True  # Оптимально

# Проверяем, является ли начальное решение оптимальным
is_optimal = check_optimality(initial_solution, c)

# Печатаем результаты
print("Initial Solution:\n", initial_solution)
print("Total Cost of Initial Solution:", total_cost)
print("Is the initial solution optimal?", is_optimal)

# Если решение не оптимально, мы обычно продолжаем с шагами оптимизации (не реализовано здесь)
