# 8.3. 오비탈과 에너지

앞에서 수소 원자의 Schrödinger 방정식의 각도 성분과 지름 성분 각각에 대한 방정식과 그 해를 살펴보았죠. 이제 둘을 합체할 차례입니다. 수소 원자의 완전한 파동 함수, 그리고 전자를 발견할 확률 분포는 어떻게 생겼는지 직접 그려서 확인해보도록 합시다.


## 지름 성분 파동 함수

8.2절에서 지름 성분의 파동 함수가 어떻게 생겼는지 대충 봤습니다. Laguerre 연관 함수가 들어있다고만 하고 넘어갔었는데요, 7장에서와 마찬가지로 처음 몇 개의 해는 직접 어떤 함수인지 보겠습니다.

|$n$|$l$|$R_{nl}(r)$ ($\rho=r/a_0$)|
|---|---|---|
| 1 | 0 | $2 \left( \frac{1}{a_0} \right)^{3/2} e^{-\rho}$ |
| 2 | 0 | $\frac{1}{2\sqrt{2}} \left( \frac{1}{a_0} \right)^{3/2} (2-\rho) e^{-\rho/2}$ |
| 2 | 1 | $\frac{1}{\sqrt{24}} \left( \frac{1}{a_0} \right)^{3/2} \rho e^{-\rho/2}$ |
| 3 | 0 | $\frac{2}{81\sqrt{3}} \left( \frac{1}{a_0} \right)^{3/2} ( 27 - 18 \rho + 2 \rho^2 ) e^{-\rho/3}$ |
| 3 | 1 | $\frac{4}{81\sqrt{6}} \left( \frac{1}{a_0} \right)^{3/2} \rho (6 - \rho) e^{-\rho/3}$ |
| 3 | 2 | $\frac{4}{81\sqrt{30}} \left( \frac{1}{a_0} \right)^{3/2} \rho^2 e^{-\rho/3}$ |

식으로만 보면 좀 복잡해보이네요. 그림으로 그려봅시다.
```python
import numpy as np
import matplotlib.pyplot as plt
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


r = np.linspace(1e-6 * a0, 25 * a0, 1000)
states = [(1, 0, "1s"), (2, 0, "2s"), (2, 1, "2p"), (3, 0, "3s"), (3, 1, "3p"), (3, 2, "3d")]

fig, ax = plt.subplots(figsize=(7, 6))
for n, l, name in states:
    ax.plot(r / a0, R(n, l, r) * a0 ** (3 / 2), lw=2, label=name)
ax.axhline(0, color="black", lw=1)
ax.set_xlim(0, 25)
ax.set_ylim(-0.15, 0.6)
ax.set_xlabel("r / a0", fontsize=10)
ax.set_ylabel("R(r)", fontsize=10)
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.show()
```
![지름 파동 함수](/assets/image-75.png)

그림으로 보니 두 가지가 눈에 띕니다.  
1) $l=0$ (s 오비탈)을 제외하고 나머지는 모두 $r=0$에서 0이 됩니다. 8.2절에서 유효 퍼텐셜에 대한 얘기를 할 때 잠깐 나왔었죠. 나머지는 모두 원점에서 퍼텐셜이 무한대로 발산하니까요.  
2) $R(r)$이 0이 되는 지점이 원점 말고 또 있습니다. 이 점에서는 전자가 발견될 확률이 0이 되고, 따라서 이 점이 바로 마디가 됩니다.

여기서 주의할 점은 규격화 상수를 계산할 때 적분을 그냥 하면 안 된다는 것입니다. 직교좌표계에서는 부피 요소가 단순하게 $dxdydz$이지만 구면좌표계에서는 부피 요소가 $r^2 \sin\theta dr d\theta d\phi$가 되거든요. 반지름이 커짐에 따라 표면적이 제곱으로 비례하기 때문입니다.


## 전자는 어디서 발견되나

마찬가지로 전자가 분포할 확률 밀도도 적분할 때 $r^2$을 곱해서 적분해야 합니다. 각도에 따른 분포까지 한번에 시각화하기에는 좀 어렵기 때문에 일단 축을 하나 잡고 그 방향으로의 확률 분포를 그리는 것이 좋습니다.

$$
P(r) = r^2 |R_{nl}(r)|^2
$$

이 값을 **지름 분포 함수(radial distribution function)**라고 부릅니다. 표면적을 곱한 것과 안 곱한 것을 한번 $R_{10}(r)$의 경우에서 비교해보죠.
```python
r = np.linspace(1e-6 * a0, 6 * a0, 1000)
R1 = R(1, 0, r)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r / a0, R1**2 / (R1**2).max(), lw=2, color="steelblue", label="|R|^2")
P1 = r**2 * R1**2
ax.plot(r / a0, P1 / P1.max(), lw=2, color="crimson", label="P(r) = r^2|R|^2")
ax.axvline(1, color="gray", ls=":", lw=1)
ax.text(1.1, 0.5, "a0", fontsize=12, color="gray")
ax.set_xlim(0, 6)
ax.set_ylim(0, 1.1)
ax.set_xlabel("r / a0")
ax.set_yticks([])
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![지름 분포 함수와 비교](/assets/image-76.png)

두 곡선의 최대점을 비교해보세요. $|R|^2$는 $r=0$에서 최대값을 가지는 반면에 $P(r)$은 $r=a_0$에서 최대값을 가지죠. 두 값이 전혀 다른 의미라는 것을 뜻합니다. $|R|^2$은 '한 점'에서의 전자 확률 밀도입니다. 하지만 전체 공간에서 보면 거리 $r$인 곳 모두에서 전자가 발견될 확률은 $P(r)$로 표현되고, 이는 거리 $r$인 모든 점들에 대한 확률 밀도의 합으로 나타나니 당연히 더 큰 값을 갖게 됩니다.  
그리고 $a_0$는 이제 Bohr가 말했던 것처럼 전자가 도는 궤도의 반지름을 뜻하는 게 아니라 전자를 발견할 확률이 최대가 되는 지점의 거리를 의미하는 것이라는 것도 알 수 있죠.

다른 경우도 같이 살펴보죠. 이번에는 최대값의 위치뿐만 아니라 위치의 평균값도 계산해보겠습니다.
```python
r = np.linspace(1e-6 * a0, 25 * a0, 1000)

fig, ax = plt.subplots(figsize=(7, 5))
for n, l, name in states:
    ax.plot(r / a0, r**2 * R(n, l, r) ** 2 * a0, lw=2, label=name)
ax.set_xlim(0, 25)
ax.set_ylim(0, 0.6)
ax.set_xlabel("r / a0")
ax.set_ylabel("P(r)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()

rw = np.linspace(1e-8 * a0, 45 * a0, 50000)

print(f"{'(n, l)':>7} {'최대 위치':>10} {'<r>':>10}")
for n, l, _ in states:
    P = rw**2 * R(n, l, rw) ** 2
    r_max = rw[np.argmax(P)] / a0
    r_mean = np.trapezoid(P * rw, rw) / a0
    print(f"{str((n, l)):>7} {r_max:9.2f} a0 {r_mean:8.2f} a0")
```
![위치의 최대값과 평균값](/assets/image-77.png)
```
 (n, l)      최대 위치        <r>
 (1, 0)      1.00 a0     1.50 a0
 (2, 0)      5.24 a0     6.00 a0
 (2, 1)      4.00 a0     5.00 a0
 (3, 0)     13.07 a0    13.50 a0
 (3, 1)     12.00 a0    12.50 a0
 (3, 2)      9.00 a0    10.50 a0
```

두 값이 꽤 차이가 나는군요. 확률 분포의 모양이 오른쪽 꼬리가 길게 늘어져 있어 평균값이 최대값보다 더 크게 나왔습니다. 그리고 $n$이 같아도 $l$이 크면 최대값과 평균값이 작아집니다. 


## 오비탈 그리기

마지막으로 수소 원자의 완전한 파동 함수를 그림으로 그려봅시다. 지름 성분 파동 함수과 각도 성분 파동 함수를 합쳐서 그리면 되죠. 물론 이 함수는 3차원 공간을 입력으로 받아서 하나의 값을 출력하기 때문에 완벽하게 그리려면 4차원 그래프를 그려야 합니다. 우리는 그렇게 할 수 없으니 일정 값을 기준으로 해서 그 값을 출력하는 3차원 공간상의 점이 어디인지 정도만 그려볼 수 있습니다. 그래도 이 그림은 꽤나 많은 정보를 우리에게 줄 수 있습니다.  
4.3절에서 3차원 상자 속 입자의 경우에 사용했던 `marching_cubes`를 이번에도 씁니다.
```python
from scipy.special import sph_harm_y
from skimage import measure


def real_Y(l, m, th, ph):
    if m == 0:
        return sph_harm_y(l, 0, th, ph).real
    elif m > 0:
        return np.sqrt(2) * (-1) ** m * sph_harm_y(l, m, th, ph).real
    else:
        return np.sqrt(2) * (-1) ** m * sph_harm_y(l, -m, th, ph).imag


def orbital_grid(n, l, m, l_box, N=100):
    g = np.linspace(-l_box, l_box, N)
    X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
    r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12 * a0
    th = np.arccos(np.clip(Z / r, -1, 1))
    ph = np.arctan2(Y, X)
    return g, R(n, l, r) * real_Y(l, m, th, ph)


def draw(ax, psi, g, frac=0.25, title=""):
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


fig = plt.figure(figsize=(13, 4))
for i, (n, l, m, l_box, name) in enumerate(
    [(1, 0, 0, 8 * a0, "1s"), (2, 0, 0, 8 * a0, "2s"), (2, 1, 1, 10 * a0, "2p_x"), (3, 2, 2, 20 * a0, "3d")], 1
):
    ax = fig.add_subplot(1, 4, i, projection="3d")
    g, psi = orbital_grid(n, l, m, l_box)
    draw(ax, psi, g, frac=0.15, title=name)

plt.tight_layout()
plt.show()
```
![수소 원자의 오비탈](/assets/image-78.png)

우리에게 정말 익숙한 그림이죠. 수소 원자의 오비탈을 그림으로 그려냈습니다. 색이 빨강과 파랑의 두 개인 이유는 파동 함수의 부호에 따라 색을 다르게 칠한 건데요, 어차피 전자를 발견할 확률은 파동 함수의 제곱에 해당해서 값은 동일하게 그려집니다. 하지만 나중에 오비탈이 겹치는 상황이 오면 그 때는 확률이 아니라 파동 함수 자체가 더해져야 해서 부호가 매우 중요해지게 됩니다.

[[TIP]]
등가면(isosurface)을 그릴 때 `frac`이라는 값을 썼는데, 이 값은 파동 함수의 최대값의 얼마의 비율을 가지는 등가면을 그릴 것이냐를 결정합니다. 위에서 말했던 것처럼 원래는 4차원 그래프를 그려야 모든 것이 보이겠지만 우리가 보는 것은 3차원이라 잘라서 봐야 하거든요. 그래서 위의 $2s$ 오비탈의 그림에서 껍질이 2개가 아니라 3개가 그려지는 것입니다. 바깥쪽 봉우리를 잘라버려서 가장자리가 2개 생겼기 때문이죠.  
보통은 전자를 발견할 확률이 90%가 되는 지점을 등가면으로 정합니다. `frac`의 값을 조절해가면서 그림을 그려보세요. 4차원까지는 아니더라도 공간에서 전자가 발견될 확률이 어떻게 분포하는지 대략적으로는 볼 수 있습니다.
[[/TIP]]


## 다음 이야기

이제 우리는 수소 원자에 대한 모든 것을 알았다고 할 수 있을까요? 사실 Schrödinger 방정식을 푸는 것만으로는 발견할 수 없는 숨은 무엇인가가 존재합니다. Schrödinger 방정식은 사실 원자핵과 전자를 그냥 전하를 띤 입자로 생각하고 푸는 것이라서 서로 간의 상호작용에 대해서는 알 수 있지만 양성자와 중성자, 전자가 어떤 구조를 가지고 있는지는 생각하지 않거든요. 전자의 구조라니 그런 게 말이 되는 건가 싶지만 놀랍게도 한 실험의 결과는 우리가 더 고려해야 할 것이 남아있다고 말해줍니다. 다음 절에서 그 얘기를 해보죠.


## 확인 문제

1. $4s$ 오비탈 ($n=4$, $l=0$)의 마디는 몇 개일까요? $4p$ 오비탈($n=4$, $l=1$)은요?
2. $1s$ 오비탈에서 전자가 $r=a_0$의 안 쪽에서 발견될 확률을 계산해보세요.
3. 오비탈 그림 코드에서 `frac`의 값을 조절해가면서 그림을 그려보세요. 본문 코드에 없는 다른 오비탈의 그림들도 한번 그려보세요.
