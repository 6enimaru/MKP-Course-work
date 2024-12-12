import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import fsolve

a = 15000
e = 0.2
mu = 398600
t = 5400


n = np.sqrt(mu / a**3)
T = 2 * np.pi / n
M = n * t

def kepler_eq(E, M, e):
    return E - e * np.sin(E) - M

E_initial_guess = M
E = fsolve(kepler_eq, E_initial_guess, args=(M, e))[0]

nu = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E / 2),
                    np.sqrt(1 - e) * np.cos(E / 2))

M_deg = np.degrees(M % (2 * np.pi))
E_deg = np.degrees(E % (2 * np.pi))
nu_deg = np.degrees(nu % (2 * np.pi))

print(f"Средняя аномалия M = {M_deg:.2f} градусов")
print(f"Эксцентриситетная аномалия E = {E_deg:.2f} градусов")
print(f"Истинная аномалия ν = {nu_deg:.2f} градусов")

times = np.linspace(0, T, 1000)
M_values = n * times

E_values = []
nu_values = []

for M_t in M_values:
    E_t = fsolve(kepler_eq, M_t, args=(M_t, e))[0]
    E_values.append(E_t)
    nu_t = 2 * np.arctan2(np.sqrt(1 + e) * np.sin(E_t / 2),
                          np.sqrt(1 - e) * np.cos(E_t / 2))
    nu_values.append(nu_t)

M_values_deg = np.degrees(np.unwrap(M_values))
E_values_deg = np.degrees(np.unwrap(E_values))
nu_values_deg = np.degrees(np.unwrap(nu_values))

#VVV графики 1 VVV
plt.figure(figsize=(12, 6))
plt.plot(times / 3600, M_values_deg, label='Средняя аномалия M')
plt.plot(times / 3600, E_values_deg, label='Эксцентриситетная аномалия E')
plt.plot(times / 3600, nu_values_deg, label='Истинная аномалия ν')
plt.xlabel('Время, часы')
plt.ylabel('Аномалия, градусы')
plt.title('Зависимость аномалий от времени')
plt.legend()
plt.grid()
plt.show()



#2часть
h = np.sqrt(mu * a * (1 - e**2))
r = a * (1 - e * np.cos(E))
V = np.sqrt(mu * (2 / r - 1 / a))
Vp = np.sqrt(2 * mu / r)
Vr = (mu / h) * e * np.sin(nu)
Vt = h / r

print(f"Радиус-вектор r = {r:.2f} км")
print(f"Модуль полной скорости V = {V:.2f} км/с")
print(f"Параболическая скорость Vp = {Vp:.2f} км/с")
print(f"Радиальная скорость Vr = {Vr:.2f} км/с")
print(f"Трансверсальная скорость Vt = {Vt:.2f} км/с")

V_values = []
Vp_values = []
Vr_values = []
Vt_values = []

for E_t, nu_t in zip(E_values, nu_values):
    r_t = a * (1 - e * np.cos(E_t))
    V_t = np.sqrt(mu * (2 / r_t - 1 / a))
    Vp_t = np.sqrt(2 * mu / r_t)
    Vr_t = (mu / h) * e * np.sin(nu_t)
    Vt_t = h / r_t

    V_values.append(V_t)
    Vp_values.append(Vp_t)
    Vr_values.append(Vr_t)
    Vt_values.append(Vt_t)

#VVV графики 2 VVV
plt.figure(figsize=(12, 6))
plt.plot(times / 3600, V_values, label='Модуль скорости V')
plt.plot(times / 3600, Vp_values, label='Параболическая скорость Vp')
plt.plot(times / 3600, Vr_values, label='Радиальная скорость Vr')
plt.plot(times / 3600, Vt_values, label='Трансверсальная скорость Vt')
plt.xlabel('Время, часы')
plt.ylabel('Скорость, км/с')
plt.title('Зависимость скоростей от времени')
plt.legend()
plt.grid()
plt.show()
