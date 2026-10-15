# 14.5. 분자 오비탈과 혼성

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/14-05.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

물 분자의 분자 오비탈 중에서 결합에 참여하지 않는 비공유 전자쌍이 들어있는 오비탈은 $1b_1$과 $3a_1$이었습니다. 둘은 서로 다른 오비탈이고 에너지도 달랐죠. 그런데 VB 이론에서 혼성 오비탈을 이야기할 때는 두 전자쌍이 동일한 것처럼 배웠습니다. 둘 중에 뭐가 맞는 걸까요?  

12.4절에서 VB 이론과 MO 이론은 서로 상보적이라고 했었죠. 같은 현상을 서로 다른 관점에서 해석하는 것이라고요. 이 말이 무슨 뜻인지 한번 확인해보도록 하겠습니다.


## 혼성 오비탈을 만드는 법

만약 두 이론이 서로 모순되는 것이 아니라면, MO 이론에서 만든 분자 오비탈로도 혼성 오비탈을 만들 수 있어야 하겠죠. 혼성 오비탈이 무엇이었는지 다시 기억해봅시다. $sp^3$ 혼성 오비탈이라는 것은 $s$ 오비탈 하나와 $p$ 오비탈 세 개가 섞인 것입니다. 이런 식으로 쓸 수 있겠네요.

$$
h_i = c_1 s + c_2 p_x + c_3 p_y + c_4 p_z
$$

총 4개의 오비탈이 섞였으니 나와야 하는 오비탈도 4개입니다. ($i=4$) 그리고 우리는 이 4개가 동등하면서도 정규직교하기를 원합니다. 일단 섞이는 비율을 생각하면 각 혼성 오비탈에서 $s$의 비율은 25%, $p$의 비율은 75%로 섞여야 하겠죠. 그러니 $|c_1|^2=\frac{1}{4}$이고, $|c_2|^2+|c_3|^2+|c_4|^2=\frac{3}{4}$가 될 겁니다. 동시에 서로 직교한다는 조건까지 고려하면 결과적으로 (회전을 무시하고) 섞는 방법은 유일하게 결정됩니다.

$$
h_1 = \frac{1}{2}(s + p_x + p_y + p_z)
$$
$$
h_2 = \frac{1}{2}(s + p_x - p_y - p_z)
$$
$$
h_3 = \frac{1}{2}(s - p_x + p_y - p_z)
$$
$$
h_4 = \frac{1}{2}(s - p_x - p_y + p_z)
$$

$p$ 오비탈은 벡터처럼 더할 수 있습니다. 예를 들어, $p_x$와 $p_y$를 반씩 더하면 $x$축과 $y$축 사이 $(1, 1)$ 방향을 향하는 $p$ 오비탈을 만들 수 있죠. 그러니까 위에서 나타낸 4개의 혼성 오비탈들은 각각 $(1, 1, 1)$, $(1, -1, -1)$, $(-1, 1, -1)$, $(-1, -1, 1)$ 방향을 향하고 있습니다. 그리고 이 네 방향은 정사면체를 이룹니다! 여기에 $s$ 오비탈이 더해지면서 부호가 같은 쪽의 로브는 커지고 반대인 쪽은 작아지죠.  

(그림 그리고 코드 지우기)
```python
import numpy as np
import matplotlib.pyplot as plt
from skimage import measure

r_H = 1.0
H1 = np.array([r_H, r_H, r_H])
H2 = np.array([r_H, -r_H, -r_H])
H3 = np.array([-r_H, r_H, -r_H])
H4 = np.array([-r_H, -r_H, r_H])

L, N = 3.5, 100
g = np.linspace(-L, L, N)
X, Y, Z = np.meshgrid(g, g, g, indexing="ij")
r = np.sqrt(X**2 + Y**2 + Z**2) + 1e-12

zeta = 1.2
s = r * np.exp(-zeta * r)
px = np.sqrt(3) * X * np.exp(-zeta * r)
py = np.sqrt(3) * Y * np.exp(-zeta * r)
pz = np.sqrt(3) * Z * np.exp(-zeta * r)

C_s = 0.4
h_1 = C_s * s + np.sqrt(1 - C_s**2) * (px + py + pz)
h_2 = C_s * s + np.sqrt(1 - C_s**2) * (px - py - pz)
h_3 = C_s * s + np.sqrt(1 - C_s**2) * (- px + py - pz)
h_4 = C_s * s + np.sqrt(1 - C_s**2) * (- px - py + pz)


def draw_molecule(ax):
    for H in (H1, H2, H3, H4):
        ax.plot([0, H[0]], [0, H[1]], [0, H[2]], color="black", lw=3)
        ax.scatter(*H, s=120, color="white", edgecolor="black", depthshade=False)
    ax.scatter(0, 0, 0, s=200, color="black", depthshade=False)


def draw_orbital(ax, psi, level, colors=("crimson", "steelblue"), alpha=0.5):
    d = g[1] - g[0]
    for lv, col in zip((level, -level), colors):
        if col is None or not psi.min() < lv < psi.max():
            continue
        v, f, _, _ = measure.marching_cubes(psi, lv, spacing=(d, d, d))
        v = v + g[0]
        ax.plot_trisurf(v[:, 0], v[:, 1], f, v[:, 2], color=col, alpha=alpha, lw=0)


def plot_setup(ax, title):
    ax.set_xlim(-L, L)
    ax.set_ylim(-L, L)
    ax.set_zlim(-L, L)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.view_init(elev=15, azim=-60)
    ax.set_title(title, fontsize=12)


level = 0.4 * np.abs(h_1).max()

fig = plt.figure(figsize=(7, 7))
panels = [
    (h_1, r"$h_1$"),
    (h_2, r"$h_2$"),
    (h_3, r"$h_3$"),
    (h_4, r"$h_4$"),
]

for k, (psi, title) in enumerate(panels, 1):
    ax = fig.add_subplot(2, 2, k, projection="3d")
    draw_orbital(ax, psi, level)
    draw_molecule(ax)
    plot_setup(ax, title)

plt.tight_layout()
plt.show()
```
![sp3 혼성 오비탈](/assets/image-105.png)

물론 모든 혼성 오비탈이 이렇게만 만들어지는 것은 아닙니다. $sp^2$나 $sp$ 혼성 오비탈의 경우에는 다른 조합으로 섞이죠. 그러니까 우리는 질문을 거꾸로 하고 있었던 겁니다. 혼성 오비탈이 되면 정사면체 구조를 갖는 것이 아니라 분자의 구조로부터 혼성 오비탈이 어떻게 만들어지는지가 결정되는 것입니다.  

물 분자의 비공유 전자쌍에 대해서도 비슷한 일을 할 수 있습니다. $1b_1$과 $3a_1$ 오비탈을 반반씩 섞어서 혼성 오비탈을 만들 수 있죠.

$$
h_+ = \frac{1}{\sqrt{2}} \left( 1b_1 + 3a_1 \right), \quad
h_- = \frac{1}{\sqrt{2}} \left( 1b_1 - 3a_1 \right)
$$

$1b_1$ 오비탈은 분자 평면에 대한 거울 반사를 하면 부호가 바뀌고, $3a_1$ 오비탈은 그대로이니 거울 반사에 대해서 $h_+$와 $h_-$는 서로 대칭입니다. 로브가 솟아있는 모습이 토끼 귀처럼 생겼다고 하여 토끼 귀 혼성(rabbit-ear hybrids)이라고도 부릅니다.


## 혼성 오비탈은 실재하는가

그런데 혼성 오비탈들이 정말로 이렇게 '동등하게' 만들어졌다면 왜 광전자 스펙트럼에서는 두 개의 에너지가 관찰되는 걸까요? Schrödinger 방정식의 형태를 다시 봅시다.  
분자 오비탈은 Hamiltonian의 고유함수여서 정해진 에너지를 가졌었죠. 원래 고유값 문제에서는 고유함수들의 선형 결합 또한 고유함수가 됩니다. 단, 고유값(에너지)이 같다는 조건에서요. 이번에 우리가 만든 혼성 오비탈은 서로 다른 에너지를 가지는 오비탈들을 섞어서 만들었습니다. 이렇게 하면 더 이상 오비탈이 Hamiltonian의 고유함수가 되지 않습니다.  

그렇다면 혼성 오비탈은 틀린 걸까요? 실제로 존재하지 않는 개념을 억지로 끼워넣어서 쓰는 걸까요? 그건 아닙니다. 혼성 오비탈은 원래의 고유함수들이 만들어내는 공간과 수학적으로 동등한 공간을 span 합니다. 즉, 원래 오비탈들로 표현할 수 있는 모든 전자 분포를 혼성 오비탈로도 완벽하게 표현할 수 있는 거죠.  

어느 쪽이 맞냐가 아니라 어떤 질문을 하느냐에 따라 편한 쪽을 쓰면 되는 겁니다. 분자의 구조나 결합을 이야기할 때는 혼성 오비탈 쪽이 편하고, 광물리학적 현상이나 전자의 들뜸과 같은 얘기를 할 때는 분자 오비탈을 쓰면 되는 거죠.


## 정리하며

이번 장에서는 실제 분자들의 분자 오비탈을 구성하는 방법에 대해 살펴봤습니다. 오비탈이 섞이는 조건 세 가지를 봤고, 분자 오비탈로부터 다양한 실험 결과들을 설명할 수 있었습니다. 그리고 마지막으로 혼성 오비탈 그림이 다른 관점에서의 해석이라는 것도 확인했습니다.  

다음 장에서는 분자 오비탈 이론을 좀 더 큰 분자에 적용하기 위해 단순화시킨 방법에 대해 알아보겠습니다. 이 방법은 특히 공액 구조를 가지는 벤젠이나 폴리메틴 염료 등의 분자들에 잘 먹힙니다.


## 확인 문제

1. $sp^2$ 혼성 오비탈은 $s$, $p_x$, $p_y$ 세 개의 오비탈을 섞어서 만듭니다. 만들어진 혼성 오비탈이 동등하고 직교하도록 계수를 정해보세요.
