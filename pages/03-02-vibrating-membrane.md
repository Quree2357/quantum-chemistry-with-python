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

여기서 만약 두 번째 항을 우변으로 넘기고 같은 논리를 쓰면 첫 번째 항과 두 번째 항도 모두 상수여야 한다는 결론을 얻습니다. 그럼 각각을 상수 $-p^2$과 $-q^2으로 두고 $p^2+q^2=k^2$이라는 관계를 만족하도록 하겠습니다.  
결과적으로 PDE 하나가 ODE 3개로 쪼개졌네요. 그리고 모두 앞에서 이미 풀어본 형태입니다.


# 양자수가 두 개라니 좋으시겠어요

경계조건을 넣어봅시다. $x$ 방향으로 하나, $y$ 방향으로 하나 있습니다. 각 방향마다 양자수가 하나씩 나오겠군요. $x$방향의 양자수를 $n$, $y$방향의 양자수를 $m$이라고 둡시다.

$$
X_n(x) = A_n \sin \left( \frac{n \pi}{a} x \right) \qquad n=1,2,3,\ldots
$$
$$
Y_m(y) = B_m \sin \left( \frac{m \pi}{b} x \right) \qquad m=1,2,3,\ldots
$$

1차원의 경우와 정확히 같은 형태지만 각 방향의 양자수가 서로 다를 수도 있으니 이번에는 아래첨자로 구분했습니다.  
각 모드의 모양은 이렇게 나옵니다. 일단 지금은 시간에 대한 함수는 잊어버리고 공간적인 모양만 보도록 합시다.  

$$
u_{nm}(x, y) = C_{nm} \sin \left( \frac{n \pi}{a} x \right) \sin \left( \frac{m \pi}{b} x \right) \qquad C_{nm} = A_n B_m
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


## 2차원의 마디는 선!

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
