# 6.1. 고전 조화진동자

분자의 결합을 용수철로 모델링하는 것을 조화진동자 모델이라고 합니다. 실제로 원자 사이에 용수철 같은 게 연결되어 있는 것은 아니지만 원자를 구성하는 원자핵과 전자 간의 전기적 힘에 의해 발생하는 진동은 용수철 모델로 잘 근사할 수 있습니다. 먼저 고전적인 조화진동자 모델을 살펴보도록 하죠. 이 모델이 왜 분자의 진동을 잘 설명할 수 있는지도 함께 확인해볼게요.


# 용수철과 에너지

용수철에 연결된 물체는 Hooke의 법칙을 따른다는 것을 배우셨을 겁니다. $F=-kx$라고 하는 그 법칙 말이죠. 여기서 $k$를 용수철 상수(spring constant) 또는 **힘 상수(force constant)**라고 불렀습니다. 단진동(simple oscillation)에 대한 Newton 운동 방정식과 그 해는 일반물리학 시간에 다뤘습니다. (잘 기억나지 않으셔도 괜찮습니다.) 질량 $m$인 물체가 힘 상수 $k$인 용수철에 매달렸을 때 물체의 운동은 다음과 같았죠.

$$
x(t)=A \cos(\omega t + \phi) \qquad \omega = \sqrt{\frac{k}{m}}
$$

$\omega$는 각진동수(angular frequency)라고 불렀고 진동수는 $\nu=\frac{\omega}{2\pi}$였습니다. 처음에 물체를 얼마나 많이 당겼는지에 따라 진폭 $A$는 다르지만 진동수는 일정했죠.  

퍼텐셜 에너지와 운동 에너지의 합으로 나타나는 전체 에너지도 보존되었습니다. 그리고 이 에너지는 진폭 $A$에 따라 결정되죠. 당연하게도 조금만 당기면 에너지가 낮고, 많이 당기면 에너지가 높습니다.

$$
E(x) = U(x) + T(x) = \frac{1}{2}kx^2 + \frac{1}{2}mv^2 = \frac{1}{2}kA^2
$$

시간에 따라서 퍼텐셜 에너지와 운동 에너지는 계속 변하지만 합은 일정하죠. 코드로 확인해볼까요?

```python
import numpy as np
import matplotlib.pyplot as plt

m, k = 1.0, 4.0
omega = np.sqrt(k / m)
A = 1.0

t = np.linspace(0, 4 * np.pi / omega, 600)
x = A * np.cos(omega * t)
v = -A * omega * np.sin(omega * t)

fig, axes = plt.subplots(1, 2, figsize=(11, 4))

axes[0].plot(t, x, color="crimson", lw=2, label="x(t)")
axes[0].axhline(0, color="black", lw=0.7)
axes[0].set_xlabel("t")
axes[0].set_ylabel("x")
axes[0].legend(fontsize=9)
axes[0].grid(alpha=0.3)

T = 0.5 * m * v**2
U = 0.5 * k * x**2

axes[1].plot(t, U, color="steelblue", lw=2, label="Potential E")
axes[1].plot(t, T, color="darkorange", lw=2, label="Kinetic E")
axes[1].plot(t, U + T, color="black", lw=2, ls="--", label="Total E")
axes[1].set_xlabel("t")
axes[1].set_ylabel("Energy")
axes[1].legend(fontsize=9)
axes[1].grid(alpha=0.3)

plt.tight_layout()
plt.show()
```
![고전 조화 진동자](/assets/image-57.png)

고전적인 조화 진동자 모델에서 중요한 점은 에너지가 연속적이라는 것입니다. 진폭 $A$가 연속적인 값을 가질 수 있으니까요. 그리고 퍼텐셜 에너지가 위치의 변화량의 제곱으로 나타난다는 점도 기억해둘 만 합니다.


## 근데 용수철이 아니잖아요

그런데 물론 처음에 얘기했던 것처럼 원자 사이에 실제로 용수철 같은 게 연결되어 있는 건 아닙니다. 그런데 왜 용수철 모델이 원자 간의 결합을 잘 근사할 수 있다고 하는 걸까요? 가장 간단한 수소 분자를 예로 들어보겠습니다.  

수소 분자는 수소 원자 두 개로 이루어져 있는 이원자 분자이고, 각각의 원자는 원자핵을 이루는 양성자 하나와 전자 하나로 구성되어 있습니다. 두 개의 수소 원자 사이의 거리에 따른 수소 분자의 퍼텐셜 에너지는 다음과 같습니다.
```python
from scipy.constants import e

De, a, re = 4.745 * e, 1.441e10, 74.16e-12
r = np.linspace(0.15e-10, 4.0e-10, 2000)
V = De * (1 - np.exp(-a * (r - re))) ** 2 - De  # 수소 분자의 Morse 근사 퍼텐셜

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(r * 1e12, V / e, color="black", lw=3)
ax.axvline(re * 1e12, color="gray", alpha=0.7, ls="--", lw=1.5)
ax.axhline(0, color='black', alpha=0.7, ls=':', lw=1.5)
ax.set_xlim(0, 400)
ax.set_ylim(-5, 2)
ax.set_xlabel("r (pm)")
ax.set_ylabel("V (eV)")
ax.grid(alpha=0.3)

plt.show()
```
![수소 분자의 퍼텐셜 에너지 곡선](/assets/image-58.png)

일반화학 시간에 많이 봤던 그래프입니다. 무한대의 거리에서 퍼텐셜 에너지가 0부터 시작해서 거리가 점점 가까워질수록 원자핵과 전자 간의 인력으로 인해 퍼텐셜 에너지가 줄어들다가 너무 가까워지면 원자핵 간의 반발력으로 인해 퍼텐셜 에너지가 다시 증가하죠. 어떻게 봐도 좌우가 비대칭이어서 용수철이랑 별로 상관이 없어보입니다.  

하지만 우리가 관심이 있는 지점은 바로 최소점입니다. 우리는 일반적으로 분자의 구조를 생각할 때 결합의 길이가 무한히 길거나 거의 0에 가까운 경우는 생각하지 않죠. 고전물리학에서와 마찬가지로 계(system)는 에너지를 최소화하는 방향으로 움직이려고 합니다. 그래서 원자들은 퍼텐셜 에너지가 가장 낮은 최소점에 위치하려고 합니다. 그래야 가장 안정해지니까요.  

퍼텐셜 에너지 곡선의 최소점 근처를 보세요. 그 주변 구간에서는 곡선이 포물선 모양으로 보입니다. 그러면 이 점 근처에서는 퍼텐셜 에너지가 위치의 변화량의 제곱에 비례해서 변할 테니 용수철 모델로 근사해도 되겠군요!  
물론 그렇게 보인다고만 하면 안 되고 실제로 그렇다는 것을 증명해야 합니다. 어떻게 해야 할까요?


## Taylor 전개 (영역 전개 아님)

임의의 함수에 대해서 어떤 한 점의 미분값들을 이용하여 원래 함수를 근사하는 것을 Taylor 전개(Taylor expansion)라고 합니다. 우리는 퍼텐셜 에너지에 관심이 있으니 거리에 대한 퍼텐셜 에너지 $V(r)$에 대한 Taylor 전개를 해봅시다. $r=r_0$인 점에서의 Taylor 전개는 다음과 같이 쓸 수 있습니다.

$$
V(r) = V(r_0) + \frac{V'(r_0)}{1!}(r-r_0) + \frac{V''(r_0)}{2!}(r-r_0)^2 + \frac{V'''(r_0)}{3!}(r-r_0)^3 + \cdots
$$

V(r_0)는 상수인데 어차피 퍼텐셜 에너지는 상수를 더하거나 빼도 상관없으니 그냥 0으로 둬도 됩니다. 또, 만약 이 점이 최소점이라면 이 점에서의 접선의 기울기, 즉 한 번 미분한 값이 0이겠죠. 다시 말해, $V'(r_0)$가 0입니다.  
그러면 남는 항이 2차항 이상이 되는데, $(r-r_0)$가 아주 작다고 가정하면 3차항 이상은 무시할 수 있습니다. 결국 최소점 근처에서는 퍼텐셜 에너지 함수를 2차 함수로 근사할 수 있는 것입니다.

$$
V(r) \approx \frac{V''(r_0)}{2}(r-r_0)^2
$$

이 식을 위의 고전 조화진동자의 경우와 비교하면 힘 상수가 다음과 같이 나타나는 것을 알 수 있습니다.

$$
k = V''(r_0) = \left.\frac{d^2}{dx^2}V(x)\right|_{x=r_0}
$$

우리는 여기서 퍼텐셜 에너지에 대한 어떠한 조건도 사용하지 않았습니다. 그냥 임의의 함수를 전개했을 뿐이죠. 그러니까 사실상 어떤 함수도 최소점 근처에서는 2차 함수로 근사할 수 있습니다! 조화진동자 모델이 물리학과 화학 곳곳에서 튀어나오는 이유가 바로 이것입니다. 수학적으로 근사하기가 좋거든요.  

실제로 겹쳐볼까요?
```python
k = 2 * De * a**2
V_harmonic = 0.5 * k * (r - re) ** 2 - De

fig, axes = plt.subplots(1, 2, figsize=(12, 5))

for ax, (xlim, ylim) in zip(axes, [((20, 300), (-5, 2)), ((60, 90), (-4.8, -4.5))]):
    ax.plot(r * 1e12, V / e, color="black", lw=3, label="actual potential")
    ax.plot(
        r * 1e12,
        V_harmonic / e,
        color="crimson",
        lw=3,
        ls="--",
        label="harmonic potential",
    )
    ax.axvline(re * 1e12, color="gray", ls=":", lw=2)
    ax.set_xlim(*xlim)
    ax.set_ylim(*ylim)
    ax.set_xlabel("r (pm)")
    ax.set_ylabel("V (eV)")
    ax.grid(alpha=0.3)

axes[0].legend(fontsize=9)
axes[0].set_title("full range")
axes[1].set_title("near minimum")
plt.tight_layout()
plt.show()
```
![퍼텐셜 에너지와 Taylor 전개](/assets/image-59.png)

왼쪽 그림에 있는 전체 곡선을 보면 실제 퍼텐셜 에너지 곡선과 포물선이 하나도 안 맞는 것처럼 보이지만 오른쪽 그림에 있는 최소점 근처에서 확대한 그림을 보면 꽤 잘 맞는 것을 볼 수 있습니다. *너무 작은 구간에서만 맞는 거 아닌가요?* 어... 사실 틀린 말은 아닌데요, 어느 정도 오차가 있는지 한번 계산해보죠.
```python
for tol in [0.01, 0.05, 0.10]:
    ok = np.abs(V_harmonic - V) < tol * De
    print(
        f"오차 {tol * 100:>2.0f}% 이내: r = {r[ok].min() * 1e12:.1f} ~ {r[ok].max() * 1e12:.1f} pm"
    )
print(f"평형 결합 길이 r_e = {re * 1e12:.1f} pm")
```
```
오차  1% 이내: r = 59.9 ~ 89.7 pm
오차  5% 이내: r = 50.4 ~ 101.5 pm
오차 10% 이내: r = 44.7 ~ 109.4 pm
평형 결합 길이 r_e = 74.2 pm
```
수소 분자의 평형 결합 길이가 74 pm 정도인데, 1% 오차 안에서 60~90 pm 범위가 커버됩니다. 결합 길이의 약 15% 정도가 변하는 거네요. 이 정도 수치면 일상적인 온도 범위에서는 벗어나지 않는 수준입니다. 물론 이보다 큰 폭으로 결합 길이가 변할 때는 조화진동사 근사가 잘 맞지 않게 되겠죠. 이 이야기는 6.4절에서 더 살펴보겠습니다.


## 다음 이야기

고전적인 조화진동자 모델을 살펴보았습니다. 그런데 이것만으로 분자의 진동을 잘 묘사하는 게 가능할까요? 사실 아직 부족한 게 맞습니다. 우리는 양자역학의 세계에 와 있으니 Schrödinger 방정식을 사용해야 하거든요. 다음 절에서는 양자역학적 조화진동자에 대해 살펴보겠습니다.


## 확인 문제

1. 힘 상수가 4배가 되면 진동수는 몇 배가 될까요? 질량이 4배가 되면요?
2. 조화 진동자의 전체 에너지가 보존되는 것을 직접 $x(t)$를 에너지 계산식에 넣어서 확인해보세요.
3. 위 코드에서 `De, a, re = 4.619 * e, 2.381, 127.5e-12`로 바꿔 HCl 분자의 경우에 퍼텐셜 에너지 곡선을 그려보세요. 수소 분자와 비교해서 곡선의 위치와 모양이 어떻게 바뀌나요?
