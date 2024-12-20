import math
import matplotlib.pyplot as plt

#дано
a = 6950
e = 0.0613
T = 5800
MU = 398600

#VVV формулы VVV
def calculate_mean_motion(mu, a):
    return math.sqrt(mu / a**3)

def calculate_mean_anomaly(n, t):
    return n * t

def calculate_eccentric_anomaly(M, e, tol=1e-6):
    E = M
    while True:
        delta_E = (M - (E - e * math.sin(E))) / (1 - e * math.cos(E))
        E += delta_E
        if abs(delta_E) < tol:
            break
    return E

def calculate_true_anomaly(E, e):
    numerator = math.sqrt(1 + e) * math.sin(E / 2)
    denominator = math.sqrt(1 - e) * math.cos(E / 2)
    nu = math.atan(numerator / denominator)
    if numerator > 0 and denominator < 0:
        nu += math.pi
    elif numerator < 0 and denominator < 0:
        nu += math.pi
    elif numerator < 0 and denominator > 0:
        nu += 2 * math.pi
    return nu

def smooth_true_anomaly(nu_values):
    smoothed_values = [nu_values[0]]
    for i in range(1, len(nu_values)):
        if nu_values[i] < smoothed_values[-1]:
            smoothed_values.append(nu_values[i] + 2 * math.pi)
        else:
            smoothed_values.append(nu_values[i])
    return smoothed_values

def calculate_radius_vector(a, e, nu):
    return a * (1 - e**2) / (1 + e * math.cos(nu))

def calculate_angular_momentum(mu, a, e):
    return math.sqrt(mu * a * (1 - e**2))

def calculate_total_velocity(mu, r, a):
    return math.sqrt(mu * (2 / r - 1 / a))

def calculate_radial_velocity(mu, h, e, nu):
    return mu / h * e * math.sin(nu)

def calculate_transverse_velocity(h, r):
    return h / r

def plot_graph(x_values, y_values, labels, title, x_label, y_label):
    plt.figure(figsize=(10, 5))
    for y, label in zip(y_values, labels):
        plt.plot(x_values, y, label=label)
    plt.title(title)
    plt.xlabel(x_label)
    plt.ylabel(y_label)
    plt.legend()
    plt.grid(True)
    plt.show()

# VVV расчеты VVV
n = calculate_mean_motion(MU, a)
M = calculate_mean_anomaly(n, T)
E = calculate_eccentric_anomaly(M, e)
nu = calculate_true_anomaly(E, e)
r = calculate_radius_vector(a, e, nu)
h = calculate_angular_momentum(MU, a, e)
V = calculate_total_velocity(MU, r, a)
Vr = calculate_radial_velocity(MU, h, e, nu)
Vt = calculate_transverse_velocity(h, r)

#VVV результаты VVV
print(f"Среднее движение (n): {n:.4f} рад/с")
print(f"Средняя аномалия (M): {M:.4f} рад")
print(f"Эксцентрическая аномалия (E): {E:.4f} рад")
print(f"Истинная аномалия (ν): {nu:.4f} рад")
print(f"Радиус-вектор (r): {r:.2f} км")
print(f"Модуль полной скорости (V): {V:.2f} км/с")
print(f"Радиальная скорость (Vr): {Vr:.2f} км/с")
print(f"Трансверсальная скорость (Vt): {Vt:.2f} км/с")

#VVV графики VVV
time_values = [i for i in range(0, int(T * 2), 10)]
M_values = [calculate_mean_anomaly(n, t) for t in time_values]
E_values = [calculate_eccentric_anomaly(M, e) for M in M_values]
nu_values = [calculate_true_anomaly(E, e) for E in E_values]
nu_values_smoothed = smooth_true_anomaly(nu_values)

r_values = [calculate_radius_vector(a, e, nu) for nu in nu_values]
V_values = [calculate_total_velocity(MU, r, a) for r in r_values]
Vr_values = [calculate_radial_velocity(MU, h, e, nu) for nu in nu_values]
Vt_values = [calculate_transverse_velocity(h, r) for r in r_values]

#VVV построение графиков VVV
plot_graph(time_values, [M_values], ["Средняя аномалия (M)"], "Средняя аномалия", "Время, сек", "Аномалия, рад")
plot_graph(time_values, [E_values], ["Эксцентрическая аномалия (E)"], "Эксцентрическая аномалия", "Время, сек", "Аномалия, рад")
plot_graph(time_values, [nu_values_smoothed], ["Истинная аномалия (ν)"], "Истинная аномалия", "Время, сек", "Аномалия, рад")
plot_graph(time_values, [r_values], ["Радиус-вектор (r)"], "Радиус-вектор", "Время, сек", "Радиус, км")
plot_graph(time_values, [V_values], ["Модуль полной скорости (V)"], "Модуль полной скорости", "Время, сек", "Скорость, км/с")
plot_graph(time_values, [Vr_values], ["Радиальная скорость (Vr)"], "Радиальная скорость", "Время, сек", "Скорость, км/с")
plot_graph(time_values, [Vt_values], ["Трансверсальная скорость (Vt)"], "Трансверсальная скорость", "Время, сек", "Скорость, км/с")
