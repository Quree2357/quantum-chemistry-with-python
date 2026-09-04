# 10.2. 침투와 오비탈 에너지

앞 절에서 원자의 에너지를 어림하기 위한 Slater 규칙에 대해 살펴보았습니다. 전자가 여러 개 있으면 핵을 가려서 유효 핵 전하가 감소할 것이라는 점에서 출발한 규칙이죠. 그런데 아직 설명하지 못하는 것이 있습니다. Slater 규칙에 따르면 $2s$ 오비탈과 $2p$ 오비탈은 같은 유효 핵 전하를 가져야 하는데 실제로는 $2s$ 오비탈이 더 안정하죠. 그 이유는 앞에서 슬쩍 힌트를 드렸는데요, 바깥쪽의 전자가 가끔 안쪽으로 파고든다고 했었던 것입니다. 그리고 그 정도가 오비탈마다 다르죠. 이번 절에서는 이것에 대해 살펴보겠습니다.


## 지름 함수 다시 쌓기 올리기

8.3절에서 만들었던 지름 함수를 다시 만들어보겠습니다. 이번에는 SciPy 패키지에서 Laguerre 함수를 직접 불러와서 임의의 $(n, l)$ 조합에 대해 쓸 수 있도록 해보죠.
```python
import numpy as np
from scipy.special import genlaguerre, factorial
from scipy.constants import physical_constants as pc

a0 = pc["Bohr radius"][0]


def R(n, l, r, Z=1):
    """수소꼴 원자의 지름 함수"""
    rho = 2 * Z * r / (n * a0)
    norm = np.sqrt((2 * Z / (n * a0))**3 * factorial(n - l - 1) / (2 * n * factorial(n + l)))
    return norm * np.exp(-rho / 2) * rho**l * genlaguerre(n - l - 1, 2 * l + 1)(rho)


def P(n, l, r, Z=1):
    """지름 분포 함수"""
    return R(n, l, r, Z) ** 2 * r**2


# 규격화 확인
r = np.linspace(1e-8 * a0, 80 * a0, 100000)
for n, l in [(2, 0), (2, 1), (3, 0), (3, 2), (4, 0)]:
    print(f"({n}, {l}): {np.trapezoid(P(n, l, r), r):.6f}")
```
```
(2, 0): 1.000000
(2, 1): 1.000000
(3, 0): 1.000000
(3, 2): 1.000000
(4, 0): 1.000000
```
8.3절에서는 Laguerre 함수를 일일이 적었는데 이제는 공식 하나로 불러올 수 있습니다. $4p$ 오비탈처럼 더 복잡한 경우도 한번에 처리할 수 있죠.


## $2s$와 $2p$가 싸우면 누가 이길까

이제 위에서 만든 코드로 리튬 원자의 $2s$ 오비탈과 $2p$ 오비탈의 지름 분포 함수를 그려봅시다. 여기서 그리는 것은 수소꼴 오비탈($Z=1$)이어서 리튬의 실제 오비탈과는 차이가 있습니다. 오비탈이 다르면 유효 핵 전하도 달라져야 하는데 $Z$ 하나로는 그걸 다 표현할 수가 없거든요. 일단 지금은 오비탈의 모양이 어떻게 달라지는가에 초점을 맞춰봅시다.
```python
import matplotlib.pyplot as plt

r = np.linspace(1e-6 * a0, 20 * a0, 2000)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(12, 5))

for n, l, name, col in [(2, 0, "2s", "crimson"), (2, 1, "2p", "steelblue")]:
    a1.plot(r / a0, P(n, l, r) * a0, lw=2, color=col, label=name)
    a2.plot(r / a0, P(n, l, r) * a0, lw=2, color=col, label=name)

a1.axvspan(0, 1, color="gray", alpha=0.2)
a1.set_xlim(0, 14)
a1.set_ylim(0, 0.2)
a1.set_xlabel("r / a0")
a1.set_ylabel("P(r)")
a1.legend(fontsize=10)
a1.grid(alpha=0.3)

a2.axvspan(0, 1, color="gray", alpha=0.2)
a2.set_xlim(0, 2)
a2.set_ylim(0, 0.08)
a2.set_xlabel("r / a0")
a2.set_ylabel("P(r)")
a2.legend(fontsize=10)
a2.grid(alpha=0.3)

plt.tight_layout()
plt.show()
```
![2s와 2p 오비탈](/assets/image-82.png)  
왼쪽 그림만 보면 $2p$ 오비탈이 더 안쪽에 있는 것처럼 보입니다. 그런데 회색으로 칠한 부분을 보세요. 핵에서 거리가 $1a_0$ 이하인 곳입니다. 리튬 원자의 경우에는 $1s$의 전자 두 개가 몰려 있는 곳이죠. 오른쪽에 확대한 그림을 보면 $2s$ 오비탈의 지름 분포 함수가 회색 부분에 봉우리를 가지고 있습니다. $2p$ 오비탈은 거의 0인데 말이죠.

실제로 계산해봅시다.
```python
r = np.linspace(1e-6 * a0, 20 * a0, 2000)
core = r < (1.0 * a0)

for n, l, name in [(2, 0, "2s"), (2, 1, "2p")]:
    inner = np.trapezoid(P(n, l, r)[core], r[core])
    r_mean = np.trapezoid(P(n, l, r) * r, r) / a0
    print(f"{name}: <r> = {r_mean:.2f} a0, r < a0일 확률: {inner * 100:.2f} %")
```
```
2s: <r> = 6.00 a0, r < a0일 확률: 3.39 %
2p: <r> = 5.00 a0, r < a0일 확률: 0.35 %
```
평균 거리만 보면 $2s$보다 $2p$가 더 가까운데, 핵에서 거리가 $1a_0$ 이하인 영역에서 전자가 발견될 확률은 거의 10배입니다. 이렇게 바깥쪽 오비탈이 안쪽 껍질을 뚫고 들어가는 현상을 **침투(penetration)**라고 합니다. 그리고 이렇게 침투가 일어나면 그동안에 전자는 가려지지 않은 핵의 전하를 그대로 느끼게 됩니다. 그러니 $2s$ 오비탈이 평균적으로 핵에 더 강하게 붙잡히고, 에너지가 더 낮게 나타나는 거죠.  

8.2절에서 수소 원자의 경우에는 $n$이 같으면 에너지도 같았는데 이제 그 축퇴가 깨졌습니다. 퍼텐셜 에너지가 $1/r$에 비례하면 우연 축퇴가 일어난다고 했었는데, 다전자 원자에서는 퍼텐셜 에너지가 다른 전자들에 의해 변형되기 때문에 더 이상 순수하게 $1/r$에 비례하지 않아 축퇴가 사라지게 됩니다.


## 3번째 껍질에서도 침투한다

3번째 껍질의 경우도 살펴봅시다. 이번에는 $3s$, $3p$, 그리고 $3d$까지 보죠.
```python
r = np.linspace(1e-6 * a0, 20 * a0, 2000)

fig, (b1, b2) = plt.subplots(1, 2, figsize=(12, 5))

for n, l, name, col in [(3, 0, "3s", "crimson"), (3, 1, "3p", "steelblue"), (3, 2, "3d", "olive")]:
    b1.plot(r / a0, P(n, l, r) * a0, lw=2, color=col, label=name)
    b2.plot(r / a0, P(n, l, r) * a0, lw=2, color=col, label=name)

b1.set_xlim(0, 20)
b1.set_ylim(0, 0.12)
b1.set_xlabel("r / a0")
b1.set_ylabel("P(r)")
b1.legend(fontsize=10)
b1.grid(alpha=0.3)

b2.set_xlim(0, 2)
b2.set_ylim(0, 0.02)
b2.set_xlabel("r / a0")
b2.set_ylabel("P(r)")
b2.legend(fontsize=10)
b2.grid(alpha=0.3)
plt.tight_layout()
plt.show()
```
![3s, 3p, 3d 오비탈](/assets/image-83.png)  
$3s$ 오비탈은 봉우리가 3개, $3p$는 2개, 그리고 $3d$ 오비탈은 1개 있습니다. 이번에도 마찬가지로 평균 거리만 보면 $s \to p \to d$ 순서로 가까워지는 것처럼 보이지만 핵 근처를 확대해보면 $3s$ 오비탈이 제일 침투를 잘 하죠.  
평균 거리와 전자 발견 확률도 다시 계산해봅시다.
```python
r = np.linspace(1e-6 * a0, 100 * a0, 2000)
core = r < (1.0 * a0)

for n, l, name in [(3, 0, "3s"), (3, 1, "3p"), (3, 2, "3d")]:
    inner = np.trapezoid(P(n, l, r)[core], r[core])
    r_mean = np.trapezoid(P(n, l, r) * r, r) / a0
    print(f"{name}: <r> = {r_mean:.2f} a0, r < a0일 확률: {inner * 100:.4f} %")
```
```
3s: <r> = 13.50 a0, r < a0일 확률: 0.9220 %
3p: <r> = 12.50 a0, r < a0일 확률: 0.1021 %
3d: <r> = 10.50 a0, r < a0일 확률: 0.0005 %
```
역시 동일한 결과가 나왔습니다. 그래서 3번째 껍질의 에너지는 $3s < 3p < 3d$순으로 커지죠. 일반화학 시간에 배웠던 그대로입니다.


## 4번째는 뭔가 다르다

마지막으로 4번째 껍질에 대해서 살펴봅시다. 그런데 여기서는 이상한 일이 벌어집니다. 혹시 포타슘(K) 원자의 전자 배치가 어떻게 됐었는지 기억하시나요? $[\text{Ar}]4s^1$이었었죠. $3d$ 오비탈이 비어 있는데도 $4s$ 오비탈에 전자가 먼저 들어갑니다. 그렇다는 얘기는 $4s$ 오비탈이 $3d$ 오비탈보다 에너지가 낮다는 건데, $n$만 보면 예상할 수 없는 일입니다. 두 오비탈의 지름 분포 함수를 그려보면 그 이유를 알 수 있습니다.
```python
r = np.linspace(1e-6 * a0, 100 * a0, 2000)

fig, ax = plt.subplots(figsize=(7, 5))

for n, l, name, col in [(4, 0, "4s", "crimson"), (3, 2, "3d", "steelblue")]:
    ax.plot(r / a0, P(n, l, r) * a0, lw=2, color=col, label=name)

ax.set_xlim(0, 30)
ax.set_ylim(0, 0.12)
ax.set_xlabel("r / a0")
ax.set_ylabel("P(r)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)
plt.show()

core = r < (1.0 * a0)

for n, l, name in [(4, 0, "4s"), (3, 2, "3d")]:
    inner = np.trapezoid(P(n, l, r)[core], r[core])
    r_mean = np.trapezoid(P(n, l, r) * r, r) / a0
    print(f"{name}: <r> = {r_mean:.2f} a0, r < a0일 확률: {inner * 100:.4f} %")
```
![](/assets/image-84.png)
```
4s: <r> = 24.00 a0, r < a0일 확률: 0.3853 %
3d: <r> = 10.50 a0, r < a0일 확률: 0.0005 %
```
$4s$ 오비탈이 $3d$ 오비탈보다 침투를 많이 해서 에너지가 낮아진 것이었네요. 전자가 발견될 확률도 거의 700배 이상 차이납니다.  

다만 지금 순서는 수소꼴 오비탈로 그렸기 때문에 실제와는 조금 다릅니다. 포타슘과 칼슘까지는 이 순서가 맞는데 그 뒤로는 두 오비탈의 에너지 차이가 크지 않아서 원자마다 순서가 바뀌는 경우가 있거든요. 실제로 전이금속의 이온화 과정을 생각해보면 분명히 중성 원자에서는 $4s$ 오비탈을 먼저 채우고 $3d$ 오비탈로 넘어갔었는데 이온화 과정에서 전자가 빠질 때는 $4s$ 전자가 먼저 빠지거든요. 일반화학 시간에 그렇게 배우긴 했지만 이유는 모르고 그냥 그렇게 된다고만 했었습니다.  
일단은 이온화 같은 건 생각하지 말고 중성 원자의 경우만 생각하도록 합시다. 어차피 우리가 다루는 대부분의 경우는 중성 원자 혹은 중성 분자니까요.


## 다음 이야기

9장에서 봤던 Pauli의 배타 원리와 더불어 이번 절에서는 오비탈의 에너지 순서가 어떻게 나타나는지 그 모양을 통해 알아보았습니다. 이제 우리는 전자가 어떻게 채워질지 예상할 수 있으니 마치 몬스터 도감을 만들듯이 주기율표를 만들어볼 수 있습니다.  각각의 성질이 어떻게 변하는지도 함께 살펴보도록 하죠.


## 확인 문제
1. $5s$ 오비탈과 $4f$ 오비탈 중 어느 쪽이 더 잘 침투할까요? 예상해본 뒤 코드로 직접 확인해보세요.
2. 위 코드에서 $Z$를 20(칼슘)으로 바꿔서 $4s$와 $3d$ 오비탈의 지름 분포 함수를 그려보세요. 핵 근처의 분포가 어떻게 달라지나요? 왜 그럴까요?
3. 수소 원자에서는 $2s$와 $2p$ 오비탈의 에너지가 같았습니다. 두 오비탈의 침투 정도는 다른데 왜 에너지가 같을까요?
