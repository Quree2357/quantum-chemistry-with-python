# 3.2. 진동하는 막 — 2차원 문제

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/03-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞 절에서 진동하는 줄에 대한 파동방정식을 풀었습니다. 이제 2차원으로 가봅시다. 줄이 아니라 막이 진동하는 겁니다.  

*Schrödinger 방정식은 언제 나오나요?* 바로 다음 절에 나올 겁니다. 조금만 기다려주세요. 지금 중요한 지점을 지나가고 있습니다. 차원을 하나 올리면 양자수가 하나 더 생기는데 이때 벌어지는 일이 나중에 중요한 열쇠가 됩니다.


## 2차원 파동방정식

2차원의 경우에는 파동방정식이 다음과 같이 됩니다. 그냥 $y$ 방향이 하나 더 추가된 거죠.

$$
\frac{\partial^2 u}{\partial x^2} + \frac{\partial^2 u}{\partial y^2} = \frac{1}{v^2} \frac{\partial^2 u}{\partial t^2}
$$

이제 함수 $u(x, y, t)$는 2차원 막 위의 각 지점이 위아래로 얼마나 움직이는지를 나타냅니다. 변수가 하나 늘었을 뿐, 의미는 똑같죠.  

가로의 길이가 $a$, 세로의 길이가 $b$인 직사각형 막을 살펴봅시다. 1차원의 경우에서와 마찬가지로 이번에도 가장자리는 고정되어 있습니다. 그러니까 다음과 같은 경계조건을 가지는 거죠.

$$
u(0, y) = u(a, y) = u(x, 0) = u(x, b) = 0 \quad \text{(모든 $t$에 대해서)}
$$


## 이번에도 변수분리!

2차원에서도 마찬가지로 변수분리법을 사용할 수 있습니다. 이번에는 변수가 3개니 분리된 함수도 3개죠.

$$
u(x, y, t)=X(x)Y(y)T(t)
$$

방정식에 대입하고 $XYT$로 양변을 나누면 이렇게 정리됩니다.

$$
\frac{1}{X} \frac{d^2 X}{dx^2} + \frac{1}{Y} \frac{d^2 Y}{dy^2} = \frac{1}{v^2 T} \frac{d^2 T}{dt^2}
$$

앞 절과 같은 논리를 쓰겠습니다. 좌변은 $x$와 $y$만의 함수이고 우변은 $t$만의 함수이니 둘 다 상수여야 합니다. 이 상수를 마찬가지로 $-k^2$으로 두도록 하죠.

그런데 이번엔 한 단계 더 나아가야 합니다. 좌변만 다시 봅시다.

$$
\frac{1}{X} \frac{d^2 X}{dx^2} + \frac{1}{Y} \frac{d^2 Y}{dy^2} = -k^2
$$

여기서 만약 두 번째 항을 우변으로 넘기고 같은 논리를 쓰면 첫 번째 항과 두 번째 항도 모두 상수여야 한다는 결론을 얻습니다. 그럼 각각을 상수 $-p^2$과 $-q^2$으로 두고 $p^2 + q^2 = k^2$이라는 관계를 만족하도록 하겠습니다.  
결과적으로 PDE 하나가 ODE 3개로 쪼개졌네요. 그리고 모두 앞에서 이미 풀어본 형태입니다.


# 양자수가 두 개라니 좋으시겠어요

경계조건을 넣어봅시다. $x$ 방향으로 하나, $y$ 방향으로 하나 있습니다. 각 방향마다 양자수가 하나씩 나오겠군요. $x$방향의 양자수를 $n$, $y$방향의 양자수를 $m$이라고 둡시다.

$$
X_n(x) = A_n \sin \left( \frac{n \pi}{a} x \right) \qquad n=1,2,3,\ldots
$$
$$
Y_m(y) = B_m \sin \left( \frac{m \pi}{b} y \right) \qquad m=1,2,3,\ldots
$$

1차원의 경우와 정확히 같은 형태지만 각 방향의 양자수가 서로 다를 수도 있으니 이번에는 아래첨자로 구분했습니다.  
각 모드의 모양은 이렇게 나옵니다. 일단 지금은 시간에 대한 함수는 잊어버리고 공간적인 모양만 보도록 합시다.  

$$
u_{nm}(x, y) = C_{nm} \sin \left( \frac{n \pi}{a} x \right) \sin \left( \frac{m \pi}{b} y \right) \qquad C_{nm} = A_n B_m
$$

음... 식으로만 봐서는 1차원 때처럼 직관적으로 모양이 느껴지지는 않습니다. 직접 그려봐야겠군요. 이번엔 2차원의 진동을 봐야 하니 3차원 그림으로 그려야 합니다.
```python
import numpy as np
import matplotlib.pyplot as plt

a, b = 2.0, 1.0
x, y = np.meshgrid(np.linspace(0, a, 100), np.linspace(0, b, 100))


def mode(n, m):
    return np.sin(n * np.pi * x / a) * np.sin(m * np.pi * y / b)


fig = plt.figure(figsize=(13, 4))

for i, (n, m) in enumerate([(1, 1), (2, 1), (2, 2), (3, 2)], 1):
    ax = fig.add_subplot(1, 4, i, projection="3d")
    ax.plot_surface(x, y, mode(n, m), cmap="coolwarm", vmin=-1, vmax=1, linewidth=0)
    ax.set_xlim(0, 2)
    ax.set_ylim(0, 1)
    ax.set_zlim(-1, 1)
    ax.set_box_aspect((1.5, 1, 1))
    ax.set_title(f"(n, m) = ({n}, {m})", fontsize=12)
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 0.5, 1])
    ax.set_zticks([-1, 0, 1])

plt.show()
```
![진동하는 막](/assets/image-35.png)

가장 간단한 모드인 (1, 1)은 가운데가 부풀어오른 모양입니다. (2, 1)에 해당하는 모드는 $x$축의 반이 갈려서 서로 반대방향으로 솟아있네요. (2, 2)는 4조각, (3, 2)는 6조각으로 나뉘어있습니다.

1차원의 진동하는 줄에서는 마디가 점으로 나타났었죠. 2차원에서는 어떻게 될까요?
```python
fig, axes = plt.subplots(1, 4, figsize=(13, 4))

for ax, (n, m) in zip(axes, [(1, 1), (2, 1), (2, 2), (3, 2)]):
    z = mode(n, m)
    ax.contourf(x, y, z, levels=20, cmap="coolwarm", vmin=-1, vmax=1)
    ax.contour(x, y, z, levels=0, colors="black", linewidth=2)
    ax.set_aspect("equal")
    ax.set_xticks([0, 1, 2])
    ax.set_yticks([0, 0.5, 1])
    ax.set_title(f"(n, m) = ({n}, {m})", fontsize=12)

plt.show()
```
![2차원 마디](/assets/image-36.png)

각 모드들을 위에서 내려다본 그림입니다. 붉은 곳이 위로, 파란 곳이 아래로 솟아있는 부분이고, 검은 선이 움직이지 않는 곳입니다.  
마디가 선이 되었네요! 이것을 마디선(nodal line)이라고 부릅니다. (2, 1)에는 세로 마디선이 하나, (2, 2)에는 가로세로 하나씩, (3, 2)에는 가로 1개, 세로 2개로 총 3개가 있네요.  

이제 3차원에서는 마디가 어떻게 나타날 지 예상이 가시죠? 네, 바로 면으로 나타납니다. 이것을 **마디면(nodal plane)**이라고 부르고요.


## 대칭성과 축퇴

자, 이제 아까 잊어버렸던 시간에 대한 함수를 다시 봅시다.

$$
\frac{1}{v^2 T} \frac{d^2 T}{dt^2} = -k^2 = -p^2 - q^2
$$

여기서 $p$와 $q$가 각각 어떻게 표현되었었는지 기억하시나요? 1차원의 경우에서와 마찬가지입니다.  

$$
p = \frac {n \pi}{a}
$$
$$
q = \frac {m \pi}{b}
$$

이제 우리는 $k$를 구할 수 있고, 시간 함수 $T(t)$도 알 수 있습니다.

$$
k_{nm} = \sqrt{p^2 + q^2} = \pi \sqrt{\frac{n^2}{a^2} + \frac{m^2}{b^2}}
$$

시간 함수는 어차피 1차원의 경우에서처럼 사인/코사인으로 나올테니 지금은 $k$값에 집중해보죠. $k$값과 파동의 진동수 $\nu$는 연관이 있었습니다.  

자, 여기서 막의 모양을 정사각형으로 만들어봅시다. 그러니까 $a=b$죠. 그러면 $k$는 다음과 같이 됩니다.

$$
k_{nm} = \frac{\pi}{a} \sqrt{n^2 + m^2}
$$

여기서 이상한 일이 벌어지는데요. (1, 2)와 (2, 1)을 한번 넣어봅시다.

$$
k_{12} = \frac{\pi}{a} \sqrt{1^2 + 2^2} = \sqrt{5} \frac{\pi}{a}
$$
$$
k_{21} = \frac{\pi}{a} \sqrt{2^2 + 1^2} = \sqrt{5} \frac{\pi}{a}
$$

어라, 두 경우의 $k$값이 같은 값을 가지는군요. 그렇다면 진동수도 같고, 에너지도 똑같다는 뜻입니다. 그런데 분명히 진동의 모양은 다릅니다.
```python
a = 1.0
x, y = np.meshgrid(np.linspace(0, a, 100), np.linspace(0, a, 100))


def mode(n, m):
    return np.sin(n * np.pi * x / a) * np.sin(m * np.pi * y / a)


fig = plt.figure(figsize=(8, 7))

for i, (n, m) in enumerate([(1, 2), (2, 1)], 1):
    ax = fig.add_subplot(2, 2, i, projection="3d")
    ax.plot_surface(x, y, mode(n, m), cmap="coolwarm", vmin=-1, vmax=1, linewidth=0)
    ax.set_title(f"(n, m) = ({n}, {m})", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([-1, 0, 1])

    ax2 = fig.add_subplot(2, 2, i + 2)
    ax2.contourf(x, y, mode(n, m), levels=20, cmap="coolwarm", vmin=-1, vmax=1)
    ax2.contour(x, y, mode(n, m), levels=0, colors="black", linewidth=2)
    ax2.set_aspect("equal")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title(f"(n, m) = ({n}, {m})", fontsize=12)

plt.tight_layout()
plt.show()
```
![2차원 축퇴](/assets/image-37.png)

이렇게 모양은 다른데 에너지(또는 진동수)가 같은 현상을 **축퇴(degeneracy)**라고 합니다. 앞으로도 자주 보게 될 개념입니다. 지금처럼 2가지의 모드가 축퇴되어 있을 때는 **이중 축퇴(doubly degenerate)**라고 합니다. two-fold degenerate라고도 하고요.

그런데 왜 이런 일이 생길까요? 막이 정사각형이기 때문입니다. 정사각형은 90도 돌려도 자기 자신과 같습니다. 그러니까 (1, 2) 모드를 90도 돌리면 (2, 1) 모드가 되는데 이 두 모드의 진동수가 다를 이유가 없습니다.  
즉, 대칭성이 축퇴를 만든다는 것입니다. 14장에서 분자의 대칭성과 군론을 배우실 때 이 이야기가 본격적으로 나올 것입니다.  

정말 정사각형에서만 일어나는 일인지 확인해봅시다. 아까처럼 직사각형의 경우에는 어떨까요?
```python
def frequency(n, m, a, b):
    return 0.5 * np.sqrt((n / a) ** 2 + (m / b) ** 2)

print("정사각형 (a, b = 1.0)")
print(f"  (1, 2) : {frequency(1, 2, 1, 1):.4f}")
print(f"  (2, 1) : {frequency(2, 1, 1, 1):.4f}")

print("직사각형 (a = 1.0, b = 2.0)")
print(f"  (1, 2) : {frequency(1, 2, 1, 2):.4f}")
print(f"  (2, 1) : {frequency(2, 1, 1, 2):.4f}")
```
```
정사각형 (a, b = 1.0)
  (1, 2) : 1.1180
  (2, 1) : 1.1180
직사각형 (a = 1.0, b = 2.0)
  (1, 2) : 0.7071
  (2, 1) : 1.0308
```
직사각형이 되니 축퇴가 깨졌습니다. 가로와 세로의 길이가 다르니 90도 회전을 하면 다른 모양이 되고 따라서 두 모드가 구별되는 것입니다. 이 현상도 나중에 보시게 될 겁니다. 원자에 자기장을 걸거나 결정 안에 존재하는 원자들의 경우에는 이렇게 대칭성이 깨지면서 에너지가 갈라지게 됩니다. 1.3절에서 "정밀하게 보면 스펙트럼 선이 갈라진다"고 했던 것도 바로 이런 이야기입니다.


## 고에너지 준위에서의 축퇴

정사각형 막에서 더 많은 모드들의 진동수를 한번 살펴보겠습니다.
```python
from collections import defaultdict

groups = defaultdict(list)

for n in range(1, 10):
    for m in range(1, 10):
        groups[n**2 + m**2].append((n, m))

for key in sorted(groups)[:10]:
    modes = groups[key]
    tag = f" ← {len(modes)}중 축퇴" if len(modes) > 1 else ""
    print(f"nu ∝ sqrt({key:2d}) = {np.sqrt(key):.4f}   {modes}{tag}")
```
```
nu ∝ sqrt( 2) = 1.4142   [(1, 1)]
nu ∝ sqrt( 5) = 2.2361   [(1, 2), (2, 1)] ← 2중 축퇴
nu ∝ sqrt( 8) = 2.8284   [(2, 2)]
nu ∝ sqrt(10) = 3.1623   [(1, 3), (3, 1)] ← 2중 축퇴
nu ∝ sqrt(13) = 3.6056   [(2, 3), (3, 2)] ← 2중 축퇴
nu ∝ sqrt(17) = 4.1231   [(1, 4), (4, 1)] ← 2중 축퇴
nu ∝ sqrt(18) = 4.2426   [(3, 3)]
nu ∝ sqrt(20) = 4.4721   [(2, 4), (4, 2)] ← 2중 축퇴
nu ∝ sqrt(25) = 5.0000   [(3, 4), (4, 3)] ← 2중 축퇴
nu ∝ sqrt(26) = 5.0990   [(1, 5), (5, 1)] ← 2중 축퇴
```
$n=m$일 때는 축퇴가 없습니다. 90도 돌려도 자기 자신이 되니까요. 나머지는 전부 이중 축퇴입니다.

그런데 더 위로 올라가면 이상한 것들이 나옵니다.
```python
for key in sorted(groups):
    modes = groups[key]
    if len(modes) > 2:
        tag = f" ← {len(modes)}중 축퇴" if len(modes) > 1 else ""
        print(f"nu ∝ sqrt({key:2d}) = {np.sqrt(key):.4f}   {modes}{tag}")
```
```
nu ∝ sqrt(50) = 7.0711   [(1, 7), (5, 5), (7, 1)] ← 3중 축퇴
nu ∝ sqrt(65) = 8.0623   [(1, 8), (4, 7), (7, 4), (8, 1)] ← 4중 축퇴
nu ∝ sqrt(85) = 9.2195   [(2, 9), (6, 7), (7, 6), (9, 2)] ← 4중 축퇴
```
고에너지 준위에서는 3중, 심지어 4중 축퇴까지 등장합니다! 이 모드들은 90도 돌린다고 겹쳐지는 것들도 아니어서 대칭성에서 오는 축퇴는 아닙니다. 이런 것을 **우연 축퇴(accidental degeneracy)**라고 부릅니다.

[[TIP]]
원형 막에서는 어떨까요? 물론 같은 방정식을 따르지만 이 경우에는 극좌표계를 사용해야 하고, Bessel 함수라는 특수함수가 튀어나와서 훨씬 해가 복잡해집니다. 여기서는 다루지 않겠지만 더 궁금하신 독자분들은 수리물리학 교재$
k_{nm} = \sqrt{p^2 + q^2} = \pi \sqrt{\frac{n^2}{a^2} + \frac{m^2}{b^2}}
$$

시간 함수는 어차피 1차원의 경우에서처럼 사인/코사인으로 나올테니 지금은 $k$값에 집중해보죠. $k$값과 파동의 진동수 $\nu$는 연관이 있었습니다.  

자, 여기서 막의 모양을 정사각형으로 만들어봅시다. 그러니까 $a=b$죠. 그러면 $k$는 다음과 같이 됩니다.

$$
k_{nm} = \frac{\pi}{a} \sqrt{n^2 + m^2}
$$

여기서 이상한 일이 벌어지는데요. (1, 2)와 (2, 1)을 한번 넣어봅시다.

$$
k_{12} = \frac{\pi}{a} \sqrt{1^2 + 2^2} = \sqrt{5} \frac{\pi}{a}
$$
$$
k_{21} = \frac{\pi}{a} \sqrt{2^2 + 1^2} = \sqrt{5} \frac{\pi}{a}
$$

어라, 두 경우의 $k$값이 같은 값을 가지는군요. 그렇다면 진동수도 같고, 에너지도 똑같다는 뜻입니다. 그런데 분명히 진동의 모양은 다릅니다.
```python
a = 1.0
x, y = np.meshgrid(np.linspace(0, a, 100), np.linspace(0, a, 100))


def mode(n, m):
    return np.sin(n * np.pi * x / a) * np.sin(m * np.pi * y / a)


fig = plt.figure(figsize=(8, 7))

for i, (n, m) in enumerate([(1, 2), (2, 1)], 1):
    ax = fig.add_subplot(2, 2, i, projection="3d")
    ax.plot_surface(x, y, mode(n, m), cmap="coolwarm", vmin=-1, vmax=1, linewidth=0)
    ax.set_title(f"(n, m) = ({n}, {m})", fontsize=12)
    ax.set_xlim(0, 1)
    ax.set_ylim(0, 1)
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([-1, 0, 1])

    ax2 = fig.add_subplot(2, 2, i + 2)
    ax2.contourf(x, y, mode(n, m), levels=20, cmap="coolwarm", vmin=-1, vmax=1)
    ax2.contour(x, y, mode(n, m), levels=0, colors="black", linewidth=2)
    ax2.set_aspect("equal")
    ax2.set_xticks([])
    ax2.set_yticks([])
    ax2.set_title(f"(n, m) = ({n}, {m})", fontsize=12)

plt.tight_layout()
plt.show()
```
![2차원 축퇴](/assets/image-37.png)

이렇게 모양은 다른데 에너지(또는 진동수)가 같은 현상을 **축퇴(degeneracy)**라고 합니다. 앞으로도 자주 보게 될 개념입니다. 지금처럼 2가지의 모드가 축퇴되어 있을 때는 **이중 축퇴(doubly degenerate)**라고 합니다. two-fold degenerate라고도 하고요.

그런데 왜 이런 일이 생길까요? 막이 정사각형이기 때문입니다. 정사각형은 90도 돌려도 자기 자신과 같습니다. 그러니까 (1, 2) 모드를 90도 돌리면 (2, 1) 모드가 되는데 이 두 모드의 진동수가 다를 이유가 없습니다.  
즉, 대칭성이 축퇴를 만든다는 것입니다. 14장에서 분자의 대칭성과 군론을 배우실 때 이 이야기가 본격적으로 나올 것입니다.  

정말 정사각형에서만 일어나는 일인지 확인해봅시다. 아까처럼 직사각형의 경우에는 어떨까요?
```python
def frequency(n, m, a, b):
    return 0.5 * np.sqrt((n / a) ** 2 + (m / b) ** 2)

print("정사각형 (a, b = 1.0)")
print(f"  (1, 2) : {frequency(1, 2, 1, 1):.4f}")
print(f"  (2, 1) : {frequency(2, 1, 1, 1):.4f}")

print("직사각형 (a = 1.0, b = 2.0)")
print(f"  (1, 2) : {frequency(1, 2, 1, 2):.4f}")
print(f"  (2, 1) : {frequency(2, 1, 1, 2):.4f}")
```
```
정사각형 (a, b = 1.0)
  (1, 2) : 1.1180
  (2, 1) : 1.1180
직사각형 (a = 1.0, b = 2.0)
  (1, 2) : 0.7071
  (2, 1) : 1.0308
```
직사각형이 되니 축퇴가 깨졌습니다. 가로와 세로의 길이가 다르니 90도 회전을 하면 다른 모양이 되고 따라서 두 모드가 구별되는 것입니다. 이 현상도 나중에 보시게 될 겁니다. 원자에 자기장을 걸거나 결정 안에 존재하는 원자들의 경우에는 이렇게 대칭성이 깨지면서 에너지가 갈라지게 됩니다. 1.3절에서 "정밀하게 보면 스펙트럼 선이 갈라진다"고 했던 것도 바로 이런 이야기입니다.


## 고에너지 준위에서의 축퇴

정사각형 막에서 더 많은 모드들의 진동수를 한번 살펴보겠습니다.
```python
from collections import defaultdict

groups = defaultdict(list)

for n in range(1, 10):
    for m in range(1, 10):
        groups[n**2 + m**2].append((n, m))

for key in sorted(groups)[:10]:
    modes = groups[key]
    tag = f" ← {len(modes)}중 축퇴" if len(modes) > 1 else ""
    print(f"nu ∝ sqrt({key:2d}) = {np.sqrt(key):.4f}   {modes}{tag}")
```
```
nu ∝ sqrt( 2) = 1.4142   [(1, 1)]
nu ∝ sqrt( 5) = 2.2361   [(1, 2), (2, 1)] ← 2중 축퇴
nu ∝ sqrt( 8) = 2.8284   [(2, 2)]
nu ∝ sqrt(10) = 3.1623   [(1, 3), (3, 1)] ← 2중 축퇴
nu ∝ sqrt(13) = 3.6056   [(2, 3), (3, 2)] ← 2중 축퇴
nu ∝ sqrt(17) = 4.1231   [(1, 4), (4, 1)] ← 2중 축퇴
nu ∝ sqrt(18) = 4.2426   [(3, 3)]
nu ∝ sqrt(20) = 4.4721   [(2, 4), (4, 2)] ← 2중 축퇴
nu ∝ sqrt(25) = 5.0000   [(3, 4), (4, 3)] ← 2중 축퇴
nu ∝ sqrt(26) = 5.0990   [(1, 5), (5, 1)] ← 2중 축퇴
```
$n=m$일 때는 축퇴가 없습니다. 90도 돌려도 자기 자신이 되니까요. 나머지는 전부 이중 축퇴입니다.

그런데 더 위로 올라가면 이상한 것들이 나옵니다.
```python
for key in sorted(groups):
    modes = groups[key]
    if len(modes) > 2:
        tag = f" ← {len(modes)}중 축퇴" if len(modes) > 1 else ""
        print(f"nu ∝ sqrt({key:2d}) = {np.sqrt(key):.4f}   {modes}{tag}")
```
```
nu ∝ sqrt(50) = 7.0711   [(1, 7), (5, 5), (7, 1)] ← 3중 축퇴
nu ∝ sqrt(65) = 8.0623   [(1, 8), (4, 7), (7, 4), (8, 1)] ← 4중 축퇴
nu ∝ sqrt(85) = 9.2195   [(2, 9), (6, 7), (7, 6), (9, 2)] ← 4중 축퇴
```
고에너지 준위에서는 3중, 심지어 4중 축퇴까지 등장합니다! 이 모드들은 90도 돌린다고 겹쳐지는 것들도 아니어서 대칭성에서 오는 축퇴는 아닙니다. 이런 것을 **우연 축퇴(accidental degeneracy)**라고 부릅니다.

[[TIP]]
원형 막에서는 어떨까요? 물론 같은 방정식을 따르지만 이 경우에는 극좌표계를 사용해야 하고, Bessel 함수라는 특수함수가 튀어나와서 훨씬 해가 복잡해집니다. 여기서는 다루지 않겠지만 더 궁금하신 독자분들은 일반적인 미분방정식 교재나 수리물리학 교재를 보시면 됩니다. 다만 원자도 구형 대칭을 가지고 있기 때문에 이 경우의 축퇴 구조는 나중에 비슷하게 써먹는 경우가 있습니다.
[[/TIP]]


## 정리하며

지금까지 1차원과 2차원에서의 고전적인 파동방정식을 풀어보았습니다. 변수분리법으로 미분방정식을 푸는 방법과 대칭성에서 오는 축퇴라는 개념도 봤죠. (물론 우연 축퇴도 있지만요.)  
이제는 *정말로* 전자에 대한 파동방정식을 알아볼 차례입니다. 그런데 전자의 경우에는 결정적으로 다른 점이 존재합니다. 과연 무엇일까요?


## 확인 문제
1. 마디선의 개수에 대한 규칙을 찾아보세요. (n, m) 모드에는 마디선이 몇 개나 있나요? (가장자리는 당연히 제외하고요.)
2. 축퇴된 두 모드를 섞은 것도 원래 파동방정식의 해가 된다는 것을 직접 대입해서 확인해보세요.
