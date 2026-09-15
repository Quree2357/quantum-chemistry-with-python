# 12.3. 화학 결합의 정체

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/12-03.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

일반화학 시간에 분자 오비탈 이론을 설명할 때 결합에 대해 얘기하면서 이 그림을 보신 기억이 있을 것입니다. 결합성 오비탈은 순수한 수소 원자의 원자 오비탈보다 에너지가 낮고 반결합성 오비탈은 에너지가 높았죠. 

![수소 분자의 MO diagram](https://chem.libretexts.org/@api/deki/files/38975/d905694ccf35a83f011afd98a63c2a6c.jpg)
[출처: https://chem.libretexts.org/Bookshelves/General_Chemistry/Map%3A_Chemistry_-_The_Central_Science_%28Brown_et_al.%29/09%3A_Molecular_Geometry_and_Bonding_Theories/9.07%3A_Molecular_Orbitals]

그런데 결합성 오비탈은 왜 에너지가 낮을까요? 그 때는 전자가 두 원자핵 사이에 몰리면서 양쪽의 원자핵과의 인력으로 안정해진다고 배웠습니다. 그럴듯한 설명입니다. 우리의 계산이 이 설명과 일치하는지 확인해보죠.


## 분자 오비탈 그리기

8.3절에서 수소 원자의 원자 오비탈을 그렸던 것과 같은 방법으로 수소 분자 이온의 분자 오비탈을 그려봅시다.

$$
\phi_{\pm} = \frac{1s_A \pm 1s_B}{\sqrt{2(1 \pm S)}}
$$

수소 원자의 $1s$ 오비탈 식을 그대로 넣고 그리면 됩니다.
```python
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure
from scipy.constants import physical_constants

a0 = physical_constants["Bohr radius"][0]


def R(n, l, r):
    rho = r / a0
    radial = {
        (1, 0): 2 * np.exp(-rho),
        (2, 0): (1 / (2 * np.sqrt(2))) * (2 - rho) * np.exp(-rho / 2),
        (2, 1): (1 / np.sqrt(24)) * rho * np.exp(-rho / 2),
        (3, 0): (2 / (81 * np.sqrt(3))) * (27 - 18 * rho + 2 * rho**2) * np.exp(-rho / 3),
        (3, 1): (4 / (81 * np.sqrt(6))) * (6 * rho - rho**2) * np.exp(-rho / 3),
        (3, 2): (4 / (81 * np.sqrt(30))) * rho**2 * np.exp(-rho / 3),
    }

    return a0 ** (-3 / 2) * radial[(n, l)]


def S(R):
    return np.exp(-R) * (1 + R + R**2 / 3)


def mo_grid(R_bohr, sign, l_box, N=100):
    R_ = R_bohr * a0
    g = np.linspace(-l_box, l_box, N)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")

    rA = np.sqrt(X**2 + Y**2 + (Z + R_ / 2) ** 2) + 1e-12 * a0
    rB = np.sqrt(X**2 + Y**2 + (Z - R_ / 2) ** 2) + 1e-12 * a0

    norm = 1 / np.sqrt(2 * (1 + sign * S(R_bohr)))
    return g, norm * (R(1, 0, rA) + sign * R(1, 0, rB))


def draw(ax, psi, g, frac, title=""):
    level = np.abs(psi).max() * frac
    d = g[1] - g[0]
    for lv, col in [(level, "crimson"), (-level, "steelblue")]:
        if not psi.min() < lv < psi.max():
            continue
        v, f, _, _ = measure.marching_cubes(psi, lv, spacing=(d,) * 3)
        v = v + g[0]
        ax.plot_trisurf(v[:, 0] / a0, v[:, 1] / a0, f, v[:, 2] / a0, color=col, alpha=0.6, lw=0)
    L = g.max() / a0
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_title(title, fontsize=12)


R_e = 2.4935

fig = plt.figure(figsize=(11, 5))

for i, (sign, name) in enumerate([(+1, "bonding"), (-1, "antibonding")], 1):
    ax = fig.add_subplot(1, 2, i, projection="3d")
    g, psi = mo_grid(R_e, sign, 5 * a0)
    draw(ax, psi, g, frac=0.2, title=name)
plt.show()
```
![수소 분자 이온의 분자 오비탈](/assets/image-92.png)

결합성 오비탈은 두 원자핵을 모두 둘러싸고 있는 하나의 덩어리처럼 보이네요. 반결합성 오비탈은 수소 원자의 $p$ 오비탈과 비슷하게 부호가 반대인 두 개의 덩어리가 나타납니다.  

단면으로 잘라보면 더 분명합니다.
```python
def mo_slice(R_bohr, sign, l_box, N=500):
    R_ = R_bohr * a0
    g = np.linspace(-l_box, l_box, N)
    X, Z = np.meshgrid(g, g)

    rA = np.sqrt(X**2 + (Z + R_ / 2) ** 2) + 1e-12 * a0
    rB = np.sqrt(X**2 + (Z - R_ / 2) ** 2) + 1e-12 * a0

    norm = 1 / np.sqrt(2 * (1 + sign * S(R_bohr)))
    return g, norm * (R(1, 0, rA) + sign * R(1, 0, rB))


fig, axes = plt.subplots(2, 2, figsize=(10, 9))

for j, (sign, name) in enumerate([(+1, "bonding"), (-1, "antibonding")]):
    g, psi = mo_slice(R_e, sign, 5 * a0)
    v = np.abs(psi).max()

    ax = axes[0, j]
    ax.pcolormesh(
        g / a0,
        g / a0,
        np.sign(psi) * np.abs(psi) ** 0.4,
        cmap="coolwarm",
        vmin=-(v**0.4),
        vmax=(v**0.4),
        shading="auto",
    )
    ax.contour(g / a0, g / a0, psi, levels=[0], colors="black", linewidths=1.5)
    ax.plot([0, 0], [-R_e / 2, R_e / 2], "o", ms=8, color="black")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{name} psi", fontsize=12)

    ax = axes[1, j]
    ax.pcolormesh(g / a0, g / a0, (psi**2)**0.4, cmap="viridis", shading="auto")
    ax.plot([0, 0], [-R_e / 2, R_e / 2], "o", ms=8, color="white")
    ax.set_aspect("equal")
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_title(f"{name} |psi|^2", fontsize=12)

plt.show()
```
![수소 분자 이온의 분자 오비탈 단면](/assets/image-93.png)

결합성 오비탈에서는 원자 오비탈 두 개가 단순하게 겹쳐진 것처럼 그려지지만 반결합성 오비탈에서는 원자핵 사이에 마디면(이 그림이 단면이라는 것을 기억하세요.)이 나타납니다. 그 면에서는 전자를 발견할 확률이 0이죠.  

4장에서 마디가 있으면 에너지가 높아진다고 했던 것 기억나시나요? 비슷하게 마디면이 있는 반결합성 오비탈의 에너지가 더 높습니다. 그리고 일반화학의 설명대로 결합성 오비탈에서는 원자핵 주변에 전자 밀도가 높게 나타나고 있네요. 이름의 유래가 바로 여기서 온 것입니다.

각각의 수소 원자 오비탈을 단순하게 더했을 때의 전자 밀도와 비교해서 그려봅시다.
```python
z = np.linspace(-4 * a0, 4 * a0, 1000)
R_ = R_e * a0
rA = np.abs(z + R_ / 2) + 1e-12 * a0
rB = np.abs(z - R_ / 2) + 1e-12 * a0
A, B = R(1, 0, rA), R(1, 0, rB)

fig, ax = plt.subplots(figsize=(7, 5))

for sign, name, col in [(+1, "bonding", "steelblue"), (-1, "antibonding", "crimson")]:
    mo = (A + sign * B) / np.sqrt(2 * (1 + sign * S(R_e)))
    ax.plot(z / a0, mo**2 * a0**3, lw=2.5, color=col, label=name)
ax.plot(z / a0, 0.5 * (A**2 + B**2) * a0**3, lw=2, ls="--", color="darkorange", label="1s_A + 1s_B")
for z in [-R_e / 2, R_e / 2]:
    ax.axvline(z, color="black", lw=1.5, ls=":")
ax.set_xlim(-3, 3)
ax.set_xlabel("z (Bohr)")
ax.set_ylabel("electron density (a.u.)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![전자 밀도 비교](/assets/image-94.png)

파란색 실선으로 그려진 결합성 오비탈과 빨간색 실선으로 그려진 반결합성 오비탈의 전자 밀도를 비교해보세요. 결합성 오비탈에서는 원자 오비탈 두 개를 단순하게 더한 것(주황색 점선)보다 가운데 부분에서 전자 밀도가 더 높습니다. 전자가 핵 근처에서 덜 발견되고 핵 사이에서 더 많이 발견된다는 뜻이죠.  


## 결합을 만드는 에너지는 무엇인가

그런데 생각해보면 전자는 음전하를 띠니까 핵 근처에 있어야 더 안정할 것 같습니다. 거리가 더 멀어진 것 같은데 그러면 퍼텐셜 에너지가 오히려 증가하지 않았을까요?

전체 에너지 말고 운동 에너지와 퍼텐셜 에너지 각각을 비교해봅시다. 다음 적분값들을 계산하면 되죠.

$$
\left\langle \hat{T} \right\rangle = \left\langle \psi_+ \middle| -\frac{1}{2}\nabla^2 \middle| \psi_+ \right\rangle \qquad \left\langle \hat{V} \right\rangle = \left\langle \psi_+ \middle| \left(-\frac{1}{r_A} -\frac{1}{r_B}+\frac{1}{R}\right) \middle| \psi_+ \right\rangle
$$

좀 복잡해보이지만 정리하고 나면 12.2절에서 계산했던 적분값들로 표현할 수 있습니다.

$$
\left\langle \hat{T} \right\rangle = \frac{1-S(R)-2K(R)}{2(1+S(R))} \qquad \left\langle \hat{V} \right\rangle = \frac{-1+J(R)+2K(R)}{1+S(R)} + \frac{1}{R}
$$

$J(R)$과 $K(R)$이라는 값은 Coulomb 적분과 공명 적분을 계산하는 과정에서 나오는 값으로 $J(R)=e^{-2R} \left(1 + \frac{1}{R}\right)-\frac{1}{R}$이고 $K(R)=-e^{-R}(1 + R)$입니다.  
분리된 원자의 경우에는 $S=0$이었고 $J$와 $K$도 0이 되므로 $\left\langle \hat{T} \right\rangle = \frac{1}{2}$이고 $\left\langle \hat{V} \right\rangle = -1$이 되겠네요. 이 값들을 기준으로 해서 에너지를 $R$에 대해 그려봅시다.
```python
def J(R):
    return np.exp(-2 * R) * (1 + 1 / R) - 1 / R


def K(R):
    return -np.exp(-R) * (1 + R)


def energies(R):
    s, j, k = S(R), J(R), K(R)
    T = (1 - s - 2 * k) / (2 + 2 * s)
    V = (-1 + j + 2 * k) / (1 + s) + 1 / R
    return T, V


R = np.linspace(1, 8, 1000)
T, V = energies(R)
i = np.argmin(T + V)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(R, T - 0.5, lw=2, color="crimson", label="delta_T")
ax.plot(R, V + 1, lw=2, color="steelblue", label="delta_V")
ax.plot(R, T + V + 0.5, lw=3, color="black", label="delta_E")
ax.axhline(0, color="gray", lw=1)
ax.axvline(R[i], color="darkorange", ls=":", lw=1.5)
ax.set_xlim(1, 8)
ax.set_ylim(-0.15, 0.2)
ax.set_xlabel("R (Bohr)")
ax.set_ylabel("Energy (hartree)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.show()

R_min = R[i]
T_min, V_min = energies(R_min)
print(f"R = {R_min:.2f}: dT = {T_min - 0.5:.4f}, dV = {V_min + 1:.4f}, dE = {T_min + V_min + 0.5:.4f}")
```
![운동 에너지와 퍼텐셜 에너지 비교](/assets/image-95.png)
```
R = 2.49: dT = -0.1173, dV = 0.0525, dE = -0.0648
```

분리된 원자($R \to \infty$)와 비교했을 때 운동 에너지는 -0.1173 hartree, 퍼텐셜 에너지는 0.0525 hartree만큼 변했습니다. 둘 중에 결합을 안정화시키는 것은 운동 에너지인 거죠. 퍼텐셜 에너지도 우리가 예상했던 것처럼 증가했습니다.  

어라, 그러면 교과서의 설명이 틀린 걸까요? 결합은 전자와 핵 사이의 정전기적 인력 때문에 안정해진다고 했으니 퍼텐셜 에너지가 낮아져야 하는데요.


## 오비탈이 수축한다

우리는 분자 오비탈을 만들 때 수소 원자의 원자 오비탈을 그대로 가져다 썼습니다. 그런데 분자를 만들고 나면 당연히 오비탈의 모양이 원래와 똑같을 리가 없죠. 매개변수 $\zeta$를 넣어서 시험 함수를 다시 만들어봅시다.

$$
\psi = c_A 1s_A (e^{-\zeta r_A}) + c_B 1s_B (e^{-\zeta r_B})
$$

$\zeta$의 값이 크면 오비탈의 크기가 작아지는 겁니다. 수소 분자가 되면 전자가 두 원자핵의 인력을 동시에 느끼게 되니 원자핵이 하나 있을 때보다 더 강하게 끌어당겨지겠죠. 그래서 크기가 수축할 것이라고 예상할 수 있습니다.

이렇게 놓고 에너지를 다시 계산해보면 다음과 같은 식을 얻습니다.

$$
\left\langle \hat{T} \right\rangle = \frac{1-S(R)-2K(R)}{2(1+S(R))} \zeta^2 \qquad \left\langle \hat{V} \right\rangle = \frac{-1+J(R)+2K(R)}{1+S(R)} \zeta + \frac{1}{R}
$$

자세한 중간 과정은 생략하겠지만 수소 분자 이온의 실제 결합 길이(약 2 Bohr)가 나타나도록 $\zeta$의 값을 최적화하면 약 1.239가 나옵니다. 에너지 그래프를 다시 그려보죠.
```python
def energies_var(zeta, R):
    s, j, k = S(zeta * R), J(zeta * R), K(zeta * R)
    T = (1 - s - 2 * k) / (2 + 2 * s) * zeta**2
    V = (-1 + j + 2 * k) / (1 + s) * zeta + 1 / R 
    return T, V


R = np.linspace(1, 100, 2000)
zeta = 1.239
T, V = energies_var(zeta, R)
i = np.argmin(T + V)

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(R, T - 0.5, lw=2, color="crimson", label="delta_T")
ax.plot(R, V + 1, lw=2, color="steelblue", label="delta_V")
ax.plot(R, T + V + 0.5, lw=3, color="black", label="delta_E")
ax.axhline(0, color="gray", lw=1)
ax.axvline(2, color="darkorange", ls=":", lw=1.5)
ax.set_xlim(1, 10)
ax.set_ylim(-0.3, 0.3)
ax.set_xlabel("R (Bohr)")
ax.set_ylabel("Energy (hartree)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.show()

R_min = R[i]
T_min, V_min = energies_var(zeta, R_min)
print(f"R = {R_min:.2f}: dT = {T_min - 0.5:.4f}, dV = {V_min + 1:.4f}, dE = {T_min + V_min + 0.5:.4f}")
```
![수축 LCAO에서의 에너지](/assets/image-96.png)
```
R = 1.99: dT = 0.0874, dV = -0.1739, dE = -0.0865
```

단순하게 더하기만 한 LCAO의 경우보다 오비탈 수축을 고려한 경우가 훨씬 더 실제 분자를 잘 설명하고 있는 겁니다. 이번에는 운동 에너지의 변화가 양수이고 퍼텐셜 에너지의 변화가 음수네요. 상자 속 입자 모델을 생각하면 오비탈이 수축한다는 것은 상자가 더 작아진다는 것이니 전자의 운동 에너지가 당연히 커질 수밖에 없습니다. 퍼텐셜 에너지는 핵과 전자의 거리가 작아져 인력이 강해짐에 따라 감소하게 되죠.

미국의 화학자 Ruedenberg는 결합이 안정해지는 이유를 단순히 퍼텐셜 에너지의 감소로만 설명할 수는 없고, 이렇게 오비탈이 수축할 수 있기 때문이라고 설명하였습니다. 보통 원자에서는 오비탈이 수축하면 퍼텐셜 에너지가 감소하는 것보다 운동 에너지의 증가가 더 커지지만 전자가 원자핵들 사이에서 '공유'되면서 파동 함수가 더 부드럽게 연결될 수 있어 운동 에너지의 증가량이 줄어든다는 것이죠.  

이렇게 전자가 공유되는 것을 **비편재화(delocalization)**라고 부릅니다. 5.3절에서 운동 에너지는 $\left\langle T \right\rangle = \frac{\left\langle p^2 \right\rangle}{2m}$로 나타남을 봤었죠. 운동량 연산자가 미분으로 나타났던 걸 생각하면 파동 함수가 급격하게 변하면 운동 에너지도 크다는 것을 알 수 있습니다. 결합성 오비탈에서는 파동 함수가 전체적으로 완만해지기 때문에 운동 에너지가 낮아지는 것입니다.

같은 논리를 반결합성 오비탈에 해보면 반대가 됩니다. 마디면이 생기면서 그 근처에서 파동 함수가 급격하게 변하니 전자의 운동 에너지가 높아지고 불안정해지는 것입니다.

[[TIP]]
일반화학 시간에 결합 차수에 대해서 배우셨을 겁니다. 그 때 결합성 오비탈에 들어 있는 전자의 개수와 반결합성 오비탈에 들어 있는 전자의 개수의 차이로 계산한다고 했었죠. 지금 하고 있는 이야기와 일맥상통합니다. 결합성 오비탈에 전자가 있으면 안정해져서 결합이 더 강하게 만들어지고 반대로 반결합성 오비탈에 전자가 있으면 불안정해져서 결합이 더 약해지는 거죠.
[[/TIP]]


## 다음 이야기

이제 수소 분자 이온에 대한 이야기는 끝입니다. 그럼 다음은 뭘까요? 당연히 전자가 두 개인 수소 분자로 가야죠. 이번에는 정확한 답을 알 수 있는 방법이 실험밖에 없습니다. 핵의 위치를 고정시킨다고 해도 삼체 문제가 되거든요. 그래도 어떻게든 최대한 정확하게 분자를 기술하는 방법을 살펴보겠습니다.


## 확인 문제
1. 결합 길이 $R$을 4 bohr, 6 bohr 등으로 바꾸면서 밀도 차이 그림을 그려보세요. 결합이 약해지면 어떻게 달라지나요?
2. 결합성 오비탈의 전자 밀도가 핵 위치에서 줄어드는 것을 규격화 상수로 설명해보세요.
3. 4.2절에서 시아닌 염료의 사슬이 길어질 수록 흡수 파장이 길어진다고 했습니다. 이 절의 비편재화와 어떻게 연결되는지 설명해보세요.
