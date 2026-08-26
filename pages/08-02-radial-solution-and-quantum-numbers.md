# 8.2. 지름 방정식과 양자수

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/08-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞 절에서 수소 원자의 Schrödinger 방정식을 살펴보았습니다. 그리고 변수분리를 통해서 지름 성분과 각도 성분의 방정식을 따로 얻었죠. 각도 성분의 파동 방정식의 해는 7장에서 봤던 구면 조화 함수로 나타났습니다. 이제 지름 성분의 파동 방정식을 풀 차례입니다. 복잡하게 생겼지만 천천히 하나씩 해결해나가면 풀 수 있습니다.


## 지름 방정식

수소 원자의 지름 방정식을 다시 한번 살펴봅시다.

$$
-\frac{\hbar^2}{2m_e r^2}\frac{d}{dr} \left( r^2 \frac{d}{dr} \right) R(r) + \left[ \frac{l(l+1)\hbar^2}{2m_e r^2} - \frac{e^2}{4\pi \varepsilon_0 r} \right] R(r) = ER(r)
$$

다시 봐도 끔찍하게 생겼군요... 일단 맨 왼쪽 항에 있는 이중 미분부터 처리해봅시다. 이렇게 치환을 해보면 어떨까요? $u(r) = rR(r)$. 한번 대입해봅시다.

$$
\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{dR}{dr}\right)=\frac{1}{r^2}\frac{d}{dr}\left(r^2\frac{d}{dr}\frac{u}{r}\right)=\frac{1}{r^2}\frac{d}{dr}\left(r^2 \left(\frac{1}{r}\frac{du}{dr}-\frac{u}{r^2}\right)\right)=\frac{1}{r^2}\frac{d}{dr}\left(r\frac{du}{dr}-u\right)
=\frac{1}{r^2}\left(r\frac{d^2u}{dr^2}\right)=\frac{1}{r}\frac{d^2u}{dr^2}
$$

치환하고 보니 복잡한 미분항이 간단한 이계 미분항으로 바뀌었습니다! 원래 지름 방정식에 넣고 전체에 $r$을 곱해서 정리하면 다음과 같이 되겠네요.

$$
-\frac{\hbar^2}{2m_e}\frac{d^2}{dr^2}u(r) + \left[ \frac{l(l+1)\hbar^2}{2m_e r^2} - \frac{e^2}{4\pi \varepsilon_0 r} \right]u(r) = Eu(r)
$$

앞에서 계속 봐왔던 1차원 Schrödinger 방정식 형태가 되었습니다. 위 식에서 $u(r)$에 곱해진 대괄호 부분을 통째로 $V_{\text{eff}}(r)$라고 하면 완전히 모양이 똑같아집니다.

$$
-\frac{\hbar^2}{2m_e}\frac{d^2}{dr^2}u(r) + V_{\text{eff}}(r)u(r) = Eu(r)
$$

경계조건도 잊지 말고 확인해야죠. $u(r)=rR(r)$이었으니까 $u(0)=0$이어야 하고, 올바른 파동 함수가 되려면 당연히 무한대의 거리에서는 $u(r) \to 0$이어야 합니다.


## 유효 퍼텐셜

위에서 대괄호 부분을 $V_{\text{eff}}(r)$라고 했었습니다. 이것을 전자가 '느끼는' 퍼텐셜이라고 하여 **유효 퍼텐셜(effective potential)**이라고 부릅니다. 잘 보면 이 항은 Coulomb 퍼텐셜에 해당하는 항과 다른 하나의 항으로 이루어져 있죠. 나머지 하나는 원심력 장벽(centrifugal barrier)이라고 부릅니다. 각운동량 때문에 생긴 항인데요, 원심력처럼 밖으로 미는(부호가 양수) 힘에 해당하는 항이라서 그런 이름이 붙었습니다. 그림으로 그리면 이렇게 생겼습니다.
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e, epsilon_0

a0 = 5.29177e-11  # Bohr 반지름
k = e**2 / (4 * np.pi * epsilon_0)

r = np.linspace(0.05, 30, 500) * a0

fig, ax = plt.subplots(figsize=(7, 5))
for l, col in zip([0, 1, 2, 3], ["black", "crimson", "darkorange", "steelblue"]):
    V_eff = (-k / r + l * (l + 1) * hbar**2 / (2 * m_e * r**2)) / e
    ax.plot(r / a0, V_eff, lw=2, color=col, label=f"l = {l}")

ax.axhline(0, color="gray", lw=1)
ax.set_xlim(0, 25)
ax.set_ylim(-15, 4)
ax.set_xlabel("r / a0")
ax.set_ylabel("V_eff (eV)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![유효 퍼텐셜](/assets/image-73.png)

$l=0$일 때만 원점에서 음의 무한대까지 내려가고 나머지는 양의 무한대로 발산하네요. 이 얘기는 전자가 $l>0$ 상태에 있으면 핵 근처에 갈 수 없다는 뜻입니다. $l=0$ 상태에 있는 전자만 핵 근처로 갈 수 있죠. 이 차이가 나중에 결정적인 역할을 하게 됩니다.


## 코드로 풀어보기

7.4절에서 만들었던 solver를 지름 방정식에다 써봅시다. 1차원이니까 그대로 풀어버릴 수 있습니다. $V$ 자리에 $V_{\text{eff}}$를 넣고 $x$ 자리에 $r$을 넣으면 되죠.
```python
def solve_1d(V, x, mass, periodic=False):
    N = len(x)
    dx = x[1] - x[0]

    D2 = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2
    if periodic:
        D2[0, -1] = D2[-1, 0] = 1 / dx**2

    H = -(hbar**2) / (2 * mass) * D2 + np.diag(V)
    E, psi = np.linalg.eigh(H)

    return E, psi / np.sqrt(dx)


N = 2000
r = np.linspace(0, 80 * a0, N + 2)[1:-1]

for l in [0, 1, 2]:
    V_eff = -k / r + l * (l + 1) * hbar**2 / (2 * m_e * r**2)
    E, u = solve_1d(V_eff, r, m_e)

    print(f"l = {l}: {np.round(E[:4] / e, 4)}")
```
```
l = 0: [-13.6003  -3.4011  -1.5117  -0.8503]
l = 1: [-3.4015 -1.5118 -0.8504 -0.5441]
l = 2: [-1.5117 -0.8504 -0.5442 -0.3701]
```

위 코드는 각 $l$에 따른 에너지 값들을 보여주고 있습니다. 7.4절에서 봤던 결과들은 $n$이나 $v$에 따라 에너지 값이 한 개씩만 나왔었는데 이번에는 조금 다르네요. 하지만 같은 값들이 계속 반복해서 나옵니다: -13.60, -3.40, -1.51. 그리고 $l$이 커질수록 하나씩 밀려서 출력되죠. 이 값들의 규칙은 놀랍게도 다음과 같습니다.

$$
E_n = -\frac{13.6}{n^2} \quad \text{(eV)}
$$

1.3절에서 Bohr의 가정으로 봤던 바로 그 식이죠. 수소 원자의 선 스펙트럼을 설명하기 위해 나왔던 그 식이요! 이번에는 Schrödinger 방정식을 풀면서 자연스럽게 나온 결과입니다.

[[TIP]]
$E_1$의 값($l=0$에서 첫 번째 값)이 -13.6003 eV로 계산되는데 실제 정확한 값은 -13.6057 eV입니다. 이 오차는 퍼텐셜이 원점에서 발산해서 그 근처에서 수치 계산의 오차가 비교적 커지기 때문입니다. 격자를 매우 촘촘하게 잡으면 오차가 점점 작아지는 것을 볼 수 있습니다. 직접 `N`의 값을 조절해보세요!
[[/TIP]]


## 양자수 삼천왕의 등장

일단 지름 성분에 대한 파동 함수의 해석해를 먼저 살펴봅시다. 6.2절에서와 마찬가지로 결과만 낼름 보죠.

$$
R_{nl}(r) = N_{nl}\left(\frac{2r}{na_0}\right)^le^{-\frac{r}{na_0}}L_{n-l-1}^{2l+1}\left(\frac{2r}{na_0}\right)
$$

조화 진동자나 구면 조화 함수 때보다는 조금 더 복잡한 꼴이네요. 체하지 않게 천천히 음미해봅시다.  
앞에 붙은 상수 $N_{nl}$은 마찬가지로 규격화 상수입니다. 그리고 거리에 대한 항이 있고, 마지막에 $L$이라는 또 이상한 함수 하나가 있습니다. 이 함수는 **Laguerre 연관 함수(associated Laguerre function)**라고 합니다. Hermite 다항식이나 Legendre 다항식처럼 이 함수도 다항식 형태죠. 중요한 점은 다항식은 항상 지수함수보다 천천히 발산한다는 점입니다. 중간에 지수함수가 껴 있으니 무한대의 거리에서는 0으로 잘 수렴하게 됩니다. 그리고 또한 지수함수 바로 앞에 있는 거리의 $l$제곱 항 때문에 $l>0$이면 $r=0$에서도 0이 됩니다. 전자가 핵에 너무 가까이 갈 수 없다는 뜻이죠.  

그리고 마지막으로 또 빼먹을 수 없는 것이 바로 양자화 조건입니다. Laguerre 함수의 첨자를 보면 $n-l-1$이 있는데, 이 값은 0 이상이어야 합니다. 따라서 이로부터 $l \leq n-1$이라는 조건이 붙죠. 위 코드 결과에서 $l$이 커질수록 에너지가 하나씩 밀리는 이유가 바로 이것입니다.

이제 지름 성분에 대한 파동 함수와 각도 성분에 대한 파동 함수를 합치면 수소 원자에 대한 완전한 파동 함수를 얻을 수 있습니다! 지금까지 등장한 양자수는 $n$, $l$, 그리고 $m$, 이렇게 총 3개죠. 7.2절에서 봤던 조건을 생각하면 세 양자수의 관계는 다음과 같습니다.

$$
n=1,2,3,\dots \quad l=0,1,\dots,n-1 \quad m=-l,-l+1,\dots,l-1,l
$$

$n$은 에너지와 직접 연관되는 양자수라서 **주양자수(principal quantum number)**라고 부릅니다. $l$은 각도 성분의 방정식에서 나오는 값이라 **각운동량 양자수(angular momentum quantum number)**라고 부르고요. $m$은 7.2절에서 봤던 것처럼 자기장에 의해 생기는 양자수라서 **자기양자수(magnetic quantum number)**라고 부릅니다.
위 관계에 따르면 $n=1$일 때는 $l=m=0$이니 $1s$ 오비탈만, $n=2$일 때는 $l=m=0$인 $2s$ 오비탈과 $l=1$에 해당하는 $2p$ 오비탈이 $m=-1, 0, 1$ 각각에 대해 $2p_x$, $2p_y$, $2p_z$의 세 가지 오비탈로 나옵니다.  
일반화학 때 배웠던 양자수와 오비탈이 왜 이렇게 존재해야 하는지를 드디어 깨달았네요! 와, 이제 무언가 보이는 듯 합니다...  


## 에너지와 축퇴

수소 원자의 전자의 에너지는 주양자수 $n$에만 의존한다는 점이 신기하네요. $l$과 $m$이 다른 값을 가져도 똑같은 에너지를 갖죠. 이렇게 축퇴된 상태는 몇 개일까요? 세어보면 각 $n$에 대해서 $l$과 $m$의 서로 다른 조합은 $n^2$개가 나옵니다. 즉, 수소 원자는 $n^2$중 축퇴를 나타내는 것입니다. 직접 그려볼까요?
```python
fig, ax = plt.subplots(figsize=(7, 5))

for n in range(1, 5):
    E = -13.6057 / n**2
    for l in range(n):
        ax.hlines(E, l - 0.35, l + 0.35, lw=2, color="steelblue")
        if l == n - 1:
            ax.text(l + 0.45, E, f"n = {n}", va="center", fontsize=9)

ax.set_xlim(-0.5, 4)
ax.set_ylim(-15, 1)
ax.set_xticks(range(4))
ax.set_xticklabels(
    ["s\n(l = 0)", "p\n(l = 1)", "d\n(l = 2)", "f\n(l = 3)"]
)
ax.axhline(0, color="gray", ls="--", lw=1)
ax.set_ylabel("E (eV)")
plt.show()
```
![수소 원자의 에너지 준위](/assets/image-74.png)

$n$이 같으면 모든 오비탈이 같은 에너지를 가집니다. $s$, $p$, $d$ 상관없이요. 근데 생각해보면 뭔가 이상합니다...

$m$에 대한 축퇴는 직관적으로 이해할 수 있습니다. 구면 대칭 때문에 $z$축을 어떻게 잡든 모양이 같으니 에너지도 같아야 하거든요. 그런데 $l$에 대한 축퇴는 이렇게 설명이 안 됩니다. 모양이 전혀 다르게 나타나거든요. $l$이 다르면 유효 퍼텐셜도 달라지니까요. 대칭성으로 설명할 수 없는 우연 축퇴가 여기서도 발생한 것입니다.  
사실 이건 수소 원자에서만 나타나는 특별한 성질인데요. 원래 $1/r$에 비례하는 Coulomb 퍼텐셜에는 숨은 대칭성이 존재합니다.(Laplace-Runge-Lenz 벡터라는 것과 관련이 있는데 이 책의 범위를 넘어서니 설명은 생략하겠습니다.) 실제로 순수한 수소가 아닌 다른 원자에서는 여러 개의 전자로 인해 퍼텐셜의 구조가 달라져 축퇴가 깨집니다.


## 다음 이야기

다음 절에서는 이번 절에서 얻은 수소 원자의 완전한 파동 함수를 직접 눈으로 확인해보겠습니다. 바로 수소 원자의 오비탈 말이죠! 처음 잡아보는 중급 몬스터인 만큼 그 보상을 마음껏 누려봅시다.


## 확인 문제

1. $n=4$인 상태는 모두 몇 개일까요? $l$과 $m$의 조합을 모두 나열해보세요.
2. 위 코드에서 `N`를 바꿔가며 에너지 계산 값이 실제 값에 얼마나 가까워지는지 직접 확인해보세요.
3. 헬륨 이온은 핵전하가 $+2e$입니다. 퍼텐셜에서 $k$를 $2k$로 바꿔서 풀어보세요. 바닥 상태의 에너지가 수소와 비교해서 몇 배인가요?
4. 원심력 장벽 항에서 $l(l+1)$ 대신 $l^2$을 쓰면 어떻게 될까요? 7.2절에서 각운동량의 크기가 $\sqrt{l(l+1)}\hbar$였던 것을 떠올려보세요.
