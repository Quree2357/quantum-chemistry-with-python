# 7.4. 프로젝트 1 — 1차원 Schrödinger 방정식 solver

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/07-04.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

지난 네 개 장에서 우리는 같은 일을 반복했습니다. 퍼텐셜을 정하고, 미분 행렬을 만들고, `eigh`를 부르고, 결과를 그렸죠. 매번 코드를 새로 쓰면서요.

이제 그럴 필요가 없게 만들 차례입니다. 임의의 퍼텐셜을 넣으면 준위와 파동함수를 돌려주는 Schrödinger solver를 만들어보겠습니다.


## 1차원 Schrödinger 방정식 solver

결과적으로 만들고 싶은 것은 이 함수입니다.

```python
E, psi = solve_1d(V, x, mass)
```

- `V` — 각 격자점에서의 퍼텐셜
- `x` — 격자점
- `mass` — 질량 (또는 회전 문제에서는 관성모멘트)
- `E, psi` — 에너지 준위와 파동함수

그리고 4장부터 7장까지 다룬 문제를 전부 이 함수 하나로 풀 수 있어야 합니다.

만드는 방법은 크게 어렵지 않습니다. 지금까지 했던 것을 정리만 해주면 되죠.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e, u, h, c


def solve_1d(V, x, mass, periodic=False):
    """
    V        : 각 격자점에서의 퍼텐셜
    x        : 격자점
    mass     : 질량, 또는 회전 문제에서는 관성모멘트
    periodic : True면 주기적 경계조건 (고리 문제)

    반환값   : 에너지 (J), 파동함수 (열 방향으로 psi[:, n])
    """
    N = len(x)
    dx = x[1] - x[0]

    D2 = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2
    if periodic:
        D2[0, -1] = D2[-1, 0] = 1 / dx**2

    H = -(hbar**2) / (2 * mass) * D2 + np.diag(V)
    E, psi = np.linalg.eigh(H)

    return E, psi / np.sqrt(dx)
```

주기적 경계조건이 두 줄로 처리됩니다. 행렬의 오른쪽 위와 왼쪽 아래 구석에 값을 넣어 첫 격자점과 마지막 격자점을 이웃으로 만드는 것이죠.

$$
D_2^{\text{고리}} = \frac{1}{h^2}\begin{pmatrix}
-2 & 1 & 0 & \cdots & \boxed{1} \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & 1 \\
\boxed{1} & 0 & 0 & 1 & -2
\end{pmatrix}
$$


## Solver 코드 검증하기

4.4절에서 강조했던 원칙입니다. 새 도구를 만들면 답을 아는 문제로 먼저 확인합니다.

```python
N = 2000

# (1) 무한 우물
L = 1e-9
x = np.linspace(0, L, N + 2)[1:-1]
E, psi = solve_1d(np.zeros(N), x, m_e)

print("=== 무한 우물 (L = 1 nm) ===")
for n in range(1, 4):
    exact = n**2 * h**2 / (8 * m_e * L**2) / e
    print(f"n = {n}: {E[n-1]/e:9.5f} eV   (해석해 {exact:.5f})")

# (2) 조화진동자
k_f, mu = 572.0, 0.504 * u
omega = np.sqrt(k_f / mu)
xh = np.linspace(-5e-11, 5e-11, N + 2)[1:-1]
E, psi = solve_1d(0.5 * k_f * xh**2, xh, mu)

print("\n=== 조화진동자 (H2) ===")
for v in range(3):
    exact = (v + 0.5) * hbar * omega / e
    print(f"v = {v}: {E[v]/e:9.5f} eV   (해석해 {exact:.5f})")

# (3) 고리
r = 1e-10
I = m_e * r**2
phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)  # 끝점 제외가 중요
E, psi = solve_1d(np.zeros(len(phi)), phi, I, periodic=True)

print("\n=== 고리 (r = 1 A) ===")
for i, m in enumerate([0, 1, -1, 2, -2]):
    exact = m**2 * hbar**2 / (2 * I) / e
    print(f"{i}번째: {E[i]/e:9.5f} eV   (m={m}, 해석해 {exact:.5f})")
```
```
=== 무한 우물 (L = 1 nm) ===
n = 1:   0.37603 eV   (해석해 0.37603)
n = 2:   1.50412 eV   (해석해 1.50412)
n = 3:   3.38427 eV   (해석해 3.38427)

=== 조화진동자 (H2) ===
v = 0:   0.27208 eV   (해석해 0.27208)
v = 1:   0.81624 eV   (해석해 0.81624)
v = 2:   1.36047 eV   (해석해 1.36039)

=== 고리 (r = 1 A) ===
0번째:  -0.00000 eV   (m=0, 해석해 0.00000)
1번째:   3.80990 eV   (m=1, 해석해 3.80998)
2번째:   3.80990 eV   (m=-1, 해석해 3.80998)
3번째:  15.23868 eV   (m=2, 해석해 15.23993)
4번째:  15.23868 eV   (m=-2, 해석해 15.23993)
```

세 경우에서 모두 계산이 제대로 되었습니다. 고리에서 $\pm m$이 짝을 이뤄 나오는 것까지 제대로 재현되네요.

[[TIP]]
고리 문제에서 `endpoint=False`가 중요합니다. $\phi = 0$과 $\phi = 2\pi$는 같은 지점이니 격자에 두 번 넣으면 안 되거든요. 상자에서 `[1:-1]`로 양 끝을 잘라냈던 것과 비슷하지만 이유는 다릅니다. 상자는 그 점의 값이 0으로 정해져 있어서 뺐고, 고리는 중복이라 뺀 것이죠.
[[/TIP]]


## 그리는 함수도 만들기

매번 같은 그림을 그리니 이것도 함수로 만들어둡시다.

```python
def plot_levels(ax, x, V, E, psi, n_levels=4, scale=0.2, xunit=1e9):
    ax.plot(x * xunit, V / e, color="black", lw=3)
    span = (E[n_levels] - E[0]) / e

    for k in range(n_levels):
        f = psi[:, k]
        if f[np.argmax(np.abs(f))] < 0:
            f = -f
        ax.hlines(E[k] / e, x[0] * xunit, x[-1] * xunit, color="gray", lw=1, ls=":")
        ax.plot(x * xunit, E[k] / e + scale * span * f / np.abs(f).max(), lw=2)

    ax.set_ylim(min(0, V.min() / e) - 0.1 * span, E[n_levels] / e)
    ax.set_xlabel("x (nm)")
    ax.set_ylabel("E (eV)")
```

이제 네 가지 문제를 한 번에 그려봅시다.

```python
N = 1500
cases = []

# 무한 우물
L = 2e-9
x1 = np.linspace(0, L, N + 2)[1:-1]
E1, p1 = solve_1d(np.zeros(N), x1, m_e)
cases.append(("infinite well", x1, np.zeros(N), E1, p1))

# 조화진동자
x2 = np.linspace(-6e-11, 6e-11, N + 2)[1:-1]
V2 = 0.5 * 572.0 * x2**2
E2, p2 = solve_1d(V2, x2, 0.504 * u)
cases.append(("harmonic", x2, V2, E2, p2))

# Morse 퍼텐셜
De, a = 37244 * 100 * h * c, 1.868e10
mu_hcl = 1.008 * 34.969 / (1.008 + 34.969) * u
x3 = np.linspace(-0.5e-10, 3e-10, N + 2)[1:-1]
V3 = De * (1 - np.exp(-a * x3)) ** 2
E3, p3 = solve_1d(V3, x3, mu_hcl)
cases.append(("Morse", x3, V3, E3, p3))

# 이중 우물
x4 = np.linspace(-1.5e-9, 1.5e-9, N + 2)[1:-1]
V4 = np.where(np.abs(np.abs(x4) - 0.6e-9) < 0.3e-9, 0.0, 3.0 * e)
E4, p4 = solve_1d(V4, x4, m_e)
cases.append(("double well", x4, V4, E4, p4))

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
for ax, (title, x, V, E, psi) in zip(axes.ravel(), cases):
    plot_levels(ax, x, V, E, psi)
    ax.set_title(title, fontsize=14)

plt.tight_layout()
plt.show()
```
![파동함수 그리기](/assets/image-69.png)

퍼텐셜 배열 하나만 바꿔서 네 문제를 풀었습니다! Solver가 제대로 작동하고 있네요. 이제 여러분은 어떤 퍼텐셜이 들어와도 물리칠 수 있는 강력한 무기를 손에 넣었습니다. 수고하셨습니다! Part II에서 뵙겠습니다.