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
D_2^{\text{고리}} = \frac{1}{dx^2}\begin{pmatrix}
-2 & 1 & 0 & \cdots & \boxed{1} \\
1 & -2 & 1 & \cdots & 0 \\
0 & 1 & -2 & \cdots & 0 \\
\vdots & \vdots & \vdots & \ddots & 1 \\
\boxed{1} & 0 & 0 & 1 & -2
\end{pmatrix}
$$


## Solver 코드 검증하기

4.4절에서 강조했던 원칙입니다. 새 도구를 만들면 답을 아는 문제로 먼저 확인합니다. 우리가 해석적으로 해를 정확히 구할 수 있었던 상자 속 입자, 조화 진동자, 고리 위 입자, 이 세 가지 경우로 살펴보죠.

```python
N = 1000

# (1) 무한 우물
L = 1e-9
x = np.linspace(0, L, N + 2)[1:-1]
E, psi = solve_1d(np.zeros(N), x, m_e)

print("=== 무한 우물 (L = 1 nm) ===")
for n in range(1, 4):
    exact = n**2 * h**2 / (8 * m_e * L**2) / e
    print(f"n = {n}: {E[n-1] / e:9.5f} eV   (해석해 {exact:.5f})")

# (2) 조화 진동자
k_f, mu = 572.0, 0.504 * u
omega = np.sqrt(k_f / mu)
xh = np.linspace(-5e-11, 5e-11, N + 2)[1:-1]
E, psi = solve_1d(0.5 * k_f * xh**2, xh, mu)

print("\n=== 조화 진동자 (H2) ===")
for v in range(3):
    exact = (v + 0.5) * hbar * omega / e
    print(f"v = {v}: {E[v] / e:9.5f} eV   (해석해 {exact:.5f})")

# (3) 고리
r = 1e-10
I = m_e * r**2
phi = np.linspace(0, 2 * np.pi, 400, endpoint=False)  # 끝점 제외가 중요
E, psi = solve_1d(np.zeros(len(phi)), phi, I, periodic=True)

print("\n=== 고리 (r = 1 A) ===")
for i, m in enumerate([0, 1, -1, 2, -2]):
    exact = m**2 * hbar**2 / (2 * I) / e
    print(f"{i}번째: {E[i] / e:9.5f} eV   (m={m}, 해석해 {exact:.5f})")
```
```
=== 무한 우물 (L = 1 nm) ===
n = 1:   0.37603 eV   (해석해 0.37603)
n = 2:   1.50412 eV   (해석해 1.50412)
n = 3:   3.38425 eV   (해석해 3.38427)

=== 조화 진동자 (H2) ===
v = 0:   0.27208 eV   (해석해 0.27208)
v = 1:   0.81623 eV   (해석해 0.81624)
v = 2:   1.36046 eV   (해석해 1.36039)

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
def plot_levels(ax, x, V, E, psi, n_levels=4, scale=0.2, xunit=1e9, xlabel="x (nm)"):
    ax.plot(x * xunit, V / e, color="black", lw=3)
    span = (E[n_levels] - E[0]) / e

    for k in range(n_levels):
        f = psi[:, k]
        if f[np.argmax(np.abs(f))] < 0:
            f = -f
        ax.hlines(E[k] / e, x[0] * xunit, x[-1] * xunit, color="gray", lw=1, ls=":")
        ax.plot(x * xunit, E[k] / e + scale * span * f / np.abs(f).max(), lw=2)

    ax.set_ylim(min(0, V.min() / e) - 0.1 * span, E[n_levels] / e)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("E (eV)")
```

이제 네 가지 문제를 한 번에 그려봅시다.

```python
N = 1000
cases = []

# (1) 무한 우물
L = 2e-9
x1 = np.linspace(0, L, N + 2)[1:-1]
V1 = np.zeros(N)
E1, p1 = solve_1d(V1, x1, m_e)
cases.append(("infinite well", x1, V1, E1, p1))

# (2) 조화 진동자
x2 = np.linspace(-6e-11, 6e-11, N + 2)[1:-1]
V2 = 0.5 * 572.0 * x2**2
E2, p2 = solve_1d(V2, x2, 0.504 * u)
cases.append(("harmonic", x2, V2, E2, p2))

# (3) Morse 퍼텐셜
De, a = 37244 * 100 * h * c, 1.868e10
mu_hcl = 1.008 * 34.969 / (1.008 + 34.969) * u
x3 = np.linspace(-0.5e-10, 3e-10, N + 2)[1:-1]
V3 = De * (1 - np.exp(-a * x3)) ** 2
E3, p3 = solve_1d(V3, x3, mu_hcl)
cases.append(("Morse", x3, V3, E3, p3))

# (4) 고리 위 입자
R = 1e-10
I = m_e * R**2
x4 = np.linspace(0, 2 * np.pi, N, endpoint=False)
V4 = np.zeros(N)
E4, p4 = solve_1d(V4, x4, I, periodic=True)

fig, axes = plt.subplots(2, 2, figsize=(12, 9))
for ax, (title, x, V, E, psi) in zip(axes.ravel(), cases):
    plot_levels(ax, x, V, E, psi)
    ax.set_title(title, fontsize=14)

plot_levels(axes[1, 1], x4, np.zeros(N), E4, p4, n_levels=5, xunit=1, xlabel="phi (rad)")
axes[1, 1].set_title("ring", fontsize=14)
axes[1, 1].set_xticks([0, np.pi, 2 * np.pi])
axes[1, 1].set_xticklabels(["0", "pi", "2pi"])

plt.tight_layout()
plt.show()
```
![파동함수 그리기](/assets/image-69.png)

퍼텐셜 배열 하나만 바꿔서 네 문제를 풀었습니다! Solver가 제대로 작동하고 있네요. 특히 오른쪽 아래의 고리 위 입자에 대한 그림을 보세요. $m=0$ 상태에 해당하는 파란색 직선을 제외하고는 각 에너지 준위별로 파동함수가 두 개씩 그려져 있습니다. 위에서 에너지를 계산했을 때와 똑같이 $\pm m$의 두 상태가 축퇴되어 있는 것이죠. 두 곡선의 위상이 90도($\pi/2$)만큼 어긋나 있는데, `eigh` 함수가 $e^{\pm i \phi}$ 대신 $\cos \phi$와 $\sin \phi$를 반환하기 때문입니다.  


## Part I을 끝내며

자, 여러분은 어떤 퍼텐셜이 들어와도 물리칠 수 있는 강력한 무기를 손에 넣었습니다. 이제 이론적으로는 양자화학 문제를 모두 풀 수 있습니다! 물론 이론적으로만 말이죠... 지금까지는 퍼텐셜이 매우 간단한 모양들이었어서 코드가 답을 금방 내놓았지만 퍼텐셜이 1차원에서 2차원만 올라가도 계산에 필요한 시간은 기하급수적으로 증가합니다. '이론상으로는' 무한대의 시간이 주어져 있다면 항상 답을 얻을 수는 있겠지만 우리는 그럴 수 없죠. 그래서 많은 사람들이 계산을 효율적으로 더 빠르게 할 수 있는 여러 가지 근사법을 개발한 것입니다. 이 내용들은 나중에 Part III나 Part IV에서 보시게 될 겁니다.  

Part II에서는 '진짜 원자'를 만나실 겁니다. 지금까지 우리가 만났던 몬스터들은 사실 하급이었던 것이죠. 우리가 얻은 무기들이 얼마나 강해졌는지, 또 앞으로 얼마나 강해질 수 있는지 알아보기 위해 중급 몬스터들이 있는 곳으로 가보겠습니다. 일단 주변 마을(카페)에 들러서 커피 한 잔과 함께 재정비를 하고 오죠! 다음 모험은 조금 더 험난해질 테니까요.


## (선택) 무기를 조금 더 강화해보기

1. 코드에서 행렬을 희소행렬로 바꿔보세요. 지금은 NumPy 패키지를 써서 $N \times N$ 행렬을 통째로 만드는데, `scipy.sparse`를 쓰면 희소 행렬로 만들 수 있어 메모리를 절약할 수 있습니다. 따라서 격자점을 좀 더 늘릴 수 있고요.
2. Solver 함수가 기댓값도 계산해주도록 만들어보세요. 위치, 운동량 등의 연산자의 기댓값을 구하는 방법은 5장에서 배웠습니다.
3. 시간에 따른 파동함수의 변화를 애니메이션으로 시각화할 수 있도록 코드를 수정해보세요. 퍼텐셜이 시간에 대해 일정하다면 Schrödinger 방정식의 해는 3.3절에서 본 것처럼 시간 독립 형태로 나타납니다. 상자 속 입자에서 $n=1$과 $n=2$를 반씩 섞은 파동함수는 시간에 따라 어떻게 움직이나요? 평형상태에서 벗어난 위치에서 시작하는 조화 진동자의 파동함수는요? 유한 장벽에서의 터널링 현상은 어떻게 나타나나요?
<details>
<summary>애니메이션화 코드 예시 보기</summary>
<div markdown="1">
    
```python
import matplotlib
from matplotlib.animation import FuncAnimation, PillowWriter


def Animate(
    x,
    V,
    E,
    psi,
    Psi0,
    T=4e-14,
    frames=360,
    xunit=1e9,
    xlabel="x (nm)",
    filename="animation.gif",
    title=""
):
    dx = x[1] - x[0]

    Psi0 = Psi0 / np.sqrt(np.trapezoid(np.abs(Psi0) ** 2, x))
    c = psi.T @ Psi0 * dx

    times = np.linspace(0, T, frames)
    prob = np.array(
        [np.abs(psi @ (c * np.exp(-1j * E * t / hbar))) ** 2 for t in times]
    )

    Vplot = V / e
    has_V = np.ptp(V) > 0
    span = np.ptp(Vplot) if has_V else 1
    base = Vplot.min()
    scale = 0.8 * span / prob.max()

    fig, ax = plt.subplots(figsize=(7, 4))
    ax.plot(x * xunit, Vplot, color="black", lw=2)
    (line,) = ax.plot([], [], color="crimson", lw=2)
    patch = [None]

    ax.set_xlim(x[0] * xunit, x[-1] * xunit)
    ax.set_ylim(base - 0.05 * span, base + 1.05 * span)
    ax.set_xlabel(xlabel)
    ax.set_ylabel("V (eV) / |psi|^2")
    ttl = ax.set_title("")

    def update(i):
        p = base + prob[i] * scale
        line.set_data(x * xunit, p)
        if patch[0]:
            patch[0].remove()
        patch[0] = ax.fill_between(x * xunit, base, p, color="crimson", alpha=0.25)
        ttl.set_text(f"{title} t = {times[i]*1e15:.2f} fs")
        return (line,)

    ani = FuncAnimation(fig, update, frames=frames, blit=False)
    ani.save(filename, writer=PillowWriter(fps=30))


# (1) 상자 속 입자에서 n=1과 n=2 상태가 섞인 상황
L, N = 2e-9, 1000
x = np.linspace(0, L, N + 2)[1:-1]
V = np.zeros(N)
E, psi = solve_1d(V, x, m_e)

Psi0 = (psi[:, 0] + psi[:, 1]) / np.sqrt(2)

Animate(x, V, E, psi, Psi0, title="particle in a box")


# (2) 조화 진동자에서 초기 파동함수가 평형상태에서 30 pm 어긋난 상황
Lh, Nh = 1e-10, 1000
k, mu = 572, 0.504 * u
omega = np.sqrt(k / mu)

xh = np.linspace(-Lh, Lh, Nh + 2)[1:-1]
Vh = 0.5 * k * xh**2
Eh, psih = solve_1d(Vh, xh, mu)

alpha = mu * omega / hbar
Psi0h = np.exp(-alpha * (xh - 3e-11) ** 2 / 2)

Animate(xh, Vh, Eh, psih, Psi0h, title="harmonic oscillator")


# (3) 유한 장벽에서의 터널링
Lt, Nt = 1e-8, 1000
xt = np.linspace(0, Lt, Nt + 2)[1:-1]
Vt = np.where(np.abs(xt - Lt / 2) < 0.2e-9, 0.6 * e, 0)
Et, psit = solve_1d(Vt, xt, m_e)

x0, sigma, k0 = 2e-9, 0.5e-9, 6e9
Psi0t = np.exp(-((xt - x0) ** 2) / (2 * sigma**2)) * np.exp(1j * k0 * xt)

Animate(xt, Vt, Et, psit, Psi0t, T=2e-14, frames=240, title="tunnelling")
```

</div>
</details>
