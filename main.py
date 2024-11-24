import numpy as np

# Определяем матрицу затрат
c = np.array([
    [40, 30, 50, 45, 70],
    [42, 28, 70, 44, 60],
    [20, 36, 62, 71, 44],
    [64, 41, 21, 37, 56],
    [34, 38, 52, 61, 48],
])

# Определяем векторы предложения и спроса
a = np.array([80, 160, 200, 70, 110])
b = np.array([80, 120, 200, 70, 110])


# Создаем начальное решение с помощью метода северо-западного угла
def north_west_corner_rule(a, b):
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
initial_solution = north_west_corner_rule(a.copy(), b.copy())

# Рассчитываем общую стоимость начального решения
total_cost = np.sum(initial_solution * c[:initial_solution.shape[0], :initial_solution.shape[1]])


# Функция для расчета потенциалов
def calculate_potentials(x, c):
    m, n = x.shape
    u = np.full(m, np.nan)  # Используем NaN для неинициализированных значений
    v = np.full(n, np.nan)

    u[0] = 0  # Устанавливаем u1 = 0

    # Шаг 1: Рассчитываем v на основе u
    for i in range(m):
        for j in range(n):
            if x[i][j] > 0 and not np.isnan(u[i]):  # Проверяем, что u[i] и x[i][j] корректны
                v[j] = c[i][j] - u[i]

    # Шаг 2: Рассчитываем u на основе v
    for j in range(n):
        for i in range(m):
            if x[i][j] > 0 and not np.isnan(v[j]):  # Проверяем, что v[j] и x[i][j] корректны
                u[i] = c[i][j] - v[j]

    return u, v


# Функция для проверки оптимальности
def check_optimality(x, c):
    m, n = x.shape
    u, v = calculate_potentials(x, c)

    for i in range(m):
        for j in range(n):
            if x[i][j] == 0:  # Проверяем только не основные переменные
                if not np.isnan(u[i]) and not np.isnan(v[j]) and (u[i] + v[j] > c[i][j]):
                    return False, u, v  # Не оптимально
    return True, u, v  # Оптимально


# Функция для поиска цикла и обновления решения
def optimize_transportation(x, c):
    is_optimal, u, v = check_optimality(x, c)

    while not is_optimal:
        # Находим неосновную переменную, которая может улучшить решение
        m, n = x.shape
        delta = np.full((m, n), np.inf)  # Инициализируем изменения
        for i in range(m):
            for j in range(n):
                if x[i][j] == 0:  # Неосновная переменная
                    if not np.isnan(u[i]) and not np.isnan(v[j]):  # Проверяем на NaN
                        delta[i][j] = u[i] + v[j] - c[i][j]

        # Находим максимальное значение delta
        max_delta = np.max(delta)
        if max_delta <= 0:
            break  # Все улучшения невыгодны

        # Находим координаты переменной, которая будет добавлена в базис
        i, j = np.unravel_index(np.argmax(delta), delta.shape)

        # Находим минимальный поток по циклу
        cycle = [(i, j)]
        visited = set(cycle)

        # Используем BFS для поиска цикла
        while True:
            found = False
            for (x_i, x_j) in cycle:
                # Проверяем соседей (верх, низ, лево, право)
                for ni, nj in [(x_i - 1, x_j), (x_i + 1, x_j), (x_i, x_j - 1), (x_i, x_j + 1)]:
                    if 0 <= ni < m and 0 <= nj < n and (ni, nj) not in visited:
                        if x[ni][nj] > 0:  # Только если это базисная переменная
                            cycle.append((ni, nj))
                            visited.add((ni, nj))
                            found = True
            if not found:
                break

        # Находим минимальный поток в цикле
        min_flow = min(x[i][j] for (i, j) in cycle if x[i][j] > 0)

        # Обновляем решение по циклу
        for (i, j) in cycle:
            if x[i][j] > 0:
                x[i][j] -= min_flow
            else:
                x[i][j] += min_flow

        # Проверяем оптимальность снова
        is_optimal, u, v = check_optimality(x, c)

    return x


# Проверяем, является ли начальное решение оптимальным
is_optimal, u, v = check_optimality(initial_solution, c)

# Печатаем результаты
print("Первоначальное решение:\n", initial_solution)
print("Минимальные затраты: ", total_cost)
print("Оптимально ли первоначальное решение?", is_optimal)

# Если решение не оптимально, продолжаем с шагами оптимизации
if not is_optimal:
    optimized_solution = optimize_transportation(initial_solution, c)
    optimized_cost = np.sum(optimized_solution * c[:optimized_solution.shape[0], :optimized_solution.shape[1]])
    print("Оптимальное решение:\n", optimized_solution)
    print("Минимальные затраты: ", optimized_cost)
else:
    print("Певоночальное решение уже оптимально.")