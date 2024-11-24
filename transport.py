import numpy as np

# Определяем матрицу затрат
cost_matrix = np.array([
    [40, 30, 50, 45, 70],
    [42, 28, 70, 44, 60],
    [20, 36, 62, 71, 44],
    [64, 41, 21, 37, 56],
    [34, 38, 52, 61, 48],
])

# Определяем векторы предложения и спроса
supply = np.array([80, 160, 200, 70, 110])
demand = np.array([80, 120, 200, 70, 110])

# Создаем начальное решение с помощью метода северо-западного угла
def north_west_corner_rule(supply, demand):
    m, n = len(supply), len(demand)
    allocation = np.zeros((m, n))

    i, j = 0, 0
    while i < m and j < n:
        if supply[i] < demand[j]:
            allocation[i][j] = supply[i]
            demand[j] -= supply[i]
            i += 1
        else:
            allocation[i][j] = demand[j]
            supply[i] -= demand[j]
            j += 1

    return allocation

# Функция для печати таблицы
def print_table(allocation, cost_matrix, supply, demand):
    m, n = allocation.shape
    print("\nТекущая таблица:")
    print("\t" + "\t".join([f"N{j + 1}" for j in range(n)]) + "\t")
    for i in range(m):
        print(f"A{i + 1}\t" + "\t".join([f"{cost_matrix[i][j]}" for j in range(n)]) + "\t")
        print(f"\t" + "\t".join([f"{int(allocation[i][j])}" for j in range(n)]) + "\t")
        print(f"\t{int(supply[i])}\t")
    print("\t" + "\t".join([str(int(demand[j])) for j in range(n)]) + "\t")
    print()

# Печатаем текущую таблицу
print("Текущая таблица перед методом северо-западного угла:")
print_table(np.zeros((len(supply), len(demand))), cost_matrix, supply, demand)

# Создаем начальное допустимое решение
initial_solution = north_west_corner_rule(supply.copy(), demand.copy())

# Печатаем начальную таблицу
print("Начальное решение после метода северо-западного угла:")
print_table(initial_solution, cost_matrix, supply, demand)

# Функция для расчета потенциалов
def calculate_potentials(allocation, cost_matrix):
    m, n = allocation.shape
    u = np.full(m, np.nan)
    v = np.full(n, np.nan)

    u[0] = 0  # Устанавливаем u1 = 0

    # Рассчитываем v на основе u
    for i in range(m):
        for j in range(n):
            if allocation[i][j] > 0 and not np.isnan(u[i]):
                v[j] = cost_matrix[i][j] - u[i]

    # Рассчитываем u на основе v
    for j in range(n):
        for i in range(m):
            if allocation[i][j] > 0 and not np.isnan(v[j]):
                u[i] = cost_matrix[i][j] - v[j]

    return u, v

# Функция для проверки оптимальности
def check_optimality(allocation, cost_matrix):
    m, n = allocation.shape
    u, v = calculate_potentials(allocation, cost_matrix)

    for i in range(m):
        for j in range(n):
            if allocation[i][j] == 0:  # Проверяем только не основные переменные
                if not np.isnan(u[i]) and not np.isnan(v[j]) and (u[i] + v[j] > cost_matrix[i][j]):
                    return False, u, v  # Не оптимально
    return True, u, v  # Оптимально

# Функция для поиска цикла и обновления решения
def optimize_transportation(allocation, cost_matrix):
    is_optimal, u, v = check_optimality(allocation, cost_matrix)

    while not is_optimal:
        # Находим неосновную переменную, которая может улучшить решение
        m, n = allocation.shape
        delta = np.full((m, n), np.inf)  # Инициализируем изменения
        for i in range(m):
            for j in range(n):
                if allocation[i][j] == 0:  # Неосновная переменная
                    if not np.isnan(u[i]) and not np.isnan(v[j]):  # Проверяем на NaN
                        delta[i][j] = u[i] + v[j] - cost_matrix[i][j]

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
                        if allocation[ni][nj] > 0:  # Только если это базисная переменная
                            cycle.append((ni, nj))
                            visited.add((ni, nj))
                            found = True
            if not found:
                break

        # Находим минимальный поток в цикле
        min_flow = min(allocation[i][j] for (i, j) in cycle if allocation[i][j] > 0)

        # Обновляем решение по циклу
        for (i, j) in cycle:
            if allocation[i][j] > 0:
                allocation[i][j] -= min_flow
            else:
                allocation[i][j] += min_flow

        # Печатаем таблицу после обновления
        print_table(allocation, cost_matrix, supply, demand)

        # Проверяем оптимальность снова
        is_optimal, u, v = check_optimality(allocation, cost_matrix)

    return allocation

# Проверяем, является ли начальное решение оптимальным
is_optimal, u, v = check_optimality(initial_solution, cost_matrix)

# Печатаем результаты
print("Начальное решение:")
print_table(initial_solution, cost_matrix, supply, demand)
print("Является ли начальное решение оптимальным?", is_optimal)

# Если решение не оптимально, продолжаем с шагами оптимизации
if not is_optimal:
    optimized_solution = optimize_transportation(initial_solution)
    print("Оптимизированное решение:")
    print_table(optimized_solution, cost_matrix, supply, demand)
else:
    print("Начальное решение уже оптимально.")

# Вычисляем минимальные затраты
def calculate_minimum_cost(allocation, cost_matrix):
    total_cost = np.sum(allocation * cost_matrix)  # Используем все элементы
    return total_cost

# Расчет и вывод минимальных затрат
minimum_cost = calculate_minimum_cost(initial_solution if is_optimal else optimized_solution, cost_matrix)
print(f"Минимальные затраты составят: {minimum_cost}")
