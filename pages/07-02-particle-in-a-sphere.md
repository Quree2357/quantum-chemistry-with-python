# 7.2. 구면 위 입자

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/07-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

이제 입자가 3차원 공간에서 움직입니다. 반지름이 $r$인 구면 위에서 입자가 돌아다닌다고 가정하는 구면 위 입자 모델에 대해 살펴보죠. 상자 속 입자에서는 1차원에서 3차원으로 넘어갈 때 크게 어려운 점이 없었는데요, 이번에는 무슨 일이 벌어질까요?


## 우리는 동그란 좌표계 밖에 못 먹어요

3차원 회전을 다루려면 구면좌표계(spherical coordinates)가 필요합니다. 원통좌표계와 거의 비슷한데요, 이번에는 $z$축 대신 $z$축과 이루는 각도 $\theta$를 사용해서 $(r, \theta, \phi)$ 형태로 좌표를 표현합니다. $r$은 원점에서부터의 거리, $\theta$는 양의 $z$축과 이루는 각도로 극각(polar angle)이라 부르고, $\phi$는 점을 $xy$-평면에 투영했을 때 양의 $x$축과 이루는 각도로 방위각(azimuthal angle)이라고 부릅니다. 직교좌표계와의 관계는 이렇게 쓸 수 있습니다.

$$
x = r\sin\theta\cos\phi, \quad y = r\sin\theta\sin\phi, \quad z = r\cos\theta \quad (0 \leq r \leq \infty, \, 0 \leq \theta \leq \pi, \, 0 \leq \phi \leq 2\pi)
$$

구면좌표계에서는 Schrödinger 방정식이 조금 복잡해지는데요, 4.3절에서 봤던 Laplacian 연산자를 먼저 표현해보겠습니다.

$$
\nabla^2 = \frac{1}{r^2}\frac{\partial}{\partial r}\left(r^2\frac{\partial}{\partial r}\right) + \frac{1}{r^2\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{r^2\sin^2\theta}\frac{\partial^2}{\partial\phi^2}
$$

어우... 제가 지금 뭘 보고 있는 거죠? 다행히 우리는 반지름이 고정되어 있으니 $r$에 대한 미분 항을 날려버릴 수 있습니다. 나머지 부분을 $1/r^2$으로 묶어서 각도 Laplcaian, 또는 Legendrian이라고 부르고 $\Lambda^2$이라고 씁니다.

$$
\Lambda^2 = \frac{1}{\sin\theta}\frac{\partial}{\partial\theta}\left(\sin\theta\frac{\partial}{\partial\theta}\right) + \frac{1}{\sin^2\theta}\frac{\partial^2}{\partial\phi^2}
$$

이 기호를 쓰면 Schrödinger 방정식을 이렇게 쓸 수 있습니다.

$$
-\frac{\hbar^2}{2I}\Lambda^2\psi(\theta,\phi) = E\psi(\theta,\phi)
$$

마찬가지로 $I$는 관성모멘트 $\mu r^2$을 나타냅니다. 그나마 보기 조금 편해졌네요.


## 또 변수분리야

두 각도 변수에 대해서 변수분리법을 써봅시다... 벌써 이게 몇 번째인지 기억도 안 나네요.

$$
\psi(\theta,\phi) = \Theta(\theta)\Phi(\phi)
$$

원래 방정식에 대입하고 정리하면 마찬가지로 각 변수만에 대한 식으로 나눌 수 있습니다.

$$
\frac{\sin\theta}{\Theta}\frac{d}{d\theta}\left(\sin\theta\frac{d\Theta}{d\theta}\right) +\frac{2IE}{\hbar^2}\sin^2\theta = m^2
$$
$$
\frac{1}{\Phi}\frac{d^2\Phi}{d\phi^2} = -m^2
$$

$\Phi$에 대한 방정식은 이미 풀어본 형태입니다. 경계조건도 똑같으니 해는 $Ae^{im\phi}$이고 $m$은 정수값이 되겠네요.  

그런데 $\Theta$에 대한 방정식은 정말 복잡해 보이네요! 계수가 상수도 아니고 어떻게 푸는지 도저히 감도 안 잡힙니다. 하지만 다행히 수학자들이 이런 형태의 미분방정식들을 이미 풀어놨습니다. 이 식에서 $\cos\theta$를 $x$로 치환하면 다음과 같은 식이 얻어지는데요.

$$
(1-x^2)\frac{d^2}{dx^2}P(x)-2x\frac{d}{dx}P(x)+\left(\frac{2IE}{\hbar^2}-\frac{m^2}{1-x^2}\right)P(x)=0
$$

$m$이 정수일 때 이 방정식을 **Legendre 방정식**이라고 합니다. 여기서 $\frac{2IE}{\hbar^2}$ 항을 일반적으로 $l(l+1)$이라고 두는데요, 두 정수 $l$과 $m$에 따라 달라지는 이 방정식의 해를 **Legendre 연관 함수(Associated Legendre polynomial)**라고 부르고 $P_l^m(x)$로 씁니다.  
저 상수 항을 $l(l+1)$이라는 모양으로 둔 이유가 궁금하실 텐데, 이렇게 정해야 Legendre 연관 함수가 $-1\leq x \leq 1$ 구간에서 무한대로 발산하지 않습니다. 그리고 또 $m$의 절대값이 $l$보다 작거나 같아야 합니다. 구간을 저렇게 정한 이유는 처음에 $\cos\theta$를 $x$로 치환했기 때문이죠.


## 구면조화함수

식이 너무 많아서 머리가 복잡해졌는데요, 간단하게 정리해봅시다. 구면좌표계에서 Schrödinger 방정식의 해는 이렇게 쓸 수 있게 됩니다. (반지름이 고정되어 있어 각도 부분만 남은 것을 보통 $Y$라고 표기합니다.)

$$
Y_l^m(\theta,\phi) = N_{lm}\,P_l^m(\cos\theta)\,e^{im\phi}
$$

앞에 붙은 $N_{lm}$은 그냥 규격화 상수입니다. 그리고 특별히 이 파동함수는 **구면조화함수(spherical harmonics)**라는 이름이 붙어있습니다. 처음 몇 개의 함수는 이렇게 생겼습니다.

$$
Y_0^0 = \sqrt{\frac{1}{4\pi}}, \quad Y_1^0 = \sqrt{\frac{3}{4\pi}}\cos\theta, \quad Y_1^{\pm1}=\mp\sqrt{\frac{3}{8\pi}}\sin\theta\,e^{\pm i\phi}
$$

구면조화함수는 SciPy 패키지에도 들어 있습니다. 직접 그림으로 그려서 어떻게 생겼는지 보도록 하죠.
```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.special import sph_harm_y


def real_Y(l, m, theta, phi):
    """실수형 구면조화함수 (그리기 편하도록)"""
    if m == 0:
        return sph_harm_y(l, 0, theta, phi).real
    elif m > 0:
        return np.sqrt(2) * (-1) ** m * sph_harm_y(l, m, theta, phi).real
    else:
        return np.sqrt(2) * (-1) ** m * sph_harm_y(l, -m, theta, phi).imag


theta = np.linspace(0, np.pi, 120)
phi = np.linspace(0, 2 * np.pi, 240)
T, P = np.meshgrid(theta, phi, indexing="ij")

states = [
    (0, 0, "Y(0,0)"),
    (1, 0, "Y(1,0)"),
    (1, 1, "Y(1,1)"),
    (2, 0, "Y(2,0)"),
    (2, 1, "Y(2,1)"),
    (2, 2, "Y(2,2)"),
]

fig = plt.figure(figsize=(13, 7))
for i, (l, m, title) in enumerate(states, 1):
    ax = fig.add_subplot(2, 3, i, projection="3d")
    Y = real_Y(l, m, T, P)
    R = np.abs(Y)
    X = R * np.sin(T) * np.cos(P)
    Yc = R * np.sin(T) * np.sin(P)
    Z = R * np.cos(T)
    colors = plt.cm.coolwarm((np.sign(Y) + 1) / 2)
    ax.plot_surface(X, Yc, Z, facecolors=colors, rstride=2, cstride=2, linewidth=0)
    lim = R.max() * 1.1
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_zlim(-lim, lim)
    ax.set_box_aspect((1, 1, 1))
    ax.set_xticks([])
    ax.set_yticks([])
    ax.set_zticks([])
    ax.set_title(title, fontsize=12)

plt.tight_layout()
plt.show()
```
![구면조화함수](/assets/image-66.png)

어디서 많이 본 것 같습니다! 바로 일반화학 시간에 봤던 s, p, d 오비탈의 모양이네요! 그러니까 구면조화함수는 오비탈의 모양을 결정하는 함수였던 것입니다.

[[TIP]]
사실 $m \neq 0$인 경우의 구면조화함수는 복소수 값을 가집니다. 하지만 같은 에너지 준위를 갖는 $\pm m$의 두 상태는 축퇴되어 있어서 선형결합을 해도 여전히 Schrödinger 방정식을 만족하죠. 그래서 실수 값을 갖도록 만든 것이 위에서 그렸던 오비탈의 식입니다.

$$
p_x \propto \frac{1}{\sqrt2}\left(Y_1^{-1} - Y_1^{1}\right), \quad p_y \propto \frac{i}{\sqrt2}\left(Y_1^{-1} + Y_1^{1}\right)
$$

그리고 당연히 구면조화함수도 직교합니다.
[[/TIP]]


## 에너지와 축퇴

에너지는 분리상수에서 바로 나옵니다.
 
$$
E_l = \frac{l(l+1)\hbar^2}{2I}, \quad l = 0, 1, 2, \ldots
$$
 
$l^2$이 아니라 $l(l+1)$입니다. 고리에서는 $m^2$이었는데 말이죠. 그리고 에너지가 $m$에 의존하지 않습니다. $l$이 같으면 $m$이 무엇이든 에너지가 같죠. $m$은 $-l$부터 $+l$까지 $2l+1$개 있으니 $l=0$은 1개, $l=1$은 3개, $l=2$는 5개입니다. 오비탈의 개수가 s가 1개, p가 3개, d가 5개인 이유가 여기에 있습니다. 
 
```python
from scipy.constants import hbar, m_e, e
 
r = 1e-10
B = hbar**2/(2*m_e*r**2)/e
 
fig, ax = plt.subplots(figsize=(7, 5))
for l in range(4):
    E = l*(l+1)*B
    g = 2*l+1
    for j in range(g):
        x0 = 0.5 + (j - (g-1)/2)*0.13
        ax.hlines(E, x0-0.05, x0+0.05, lw=3, color="navy")
    ax.text(1.05, E, f"l = {l}   (x{g})", va="center", fontsize=9)
 
ax.set_xlim(-0.1, 1.6); ax.set_ylim(-3, 52)
ax.set_xticks([]); ax.set_ylabel("E (eV)")
plt.show()
 
for l in range(4):
    print(f"l={l}: E = {l*(l+1)*B:8.4f} eV,  m = {list(range(-l, l+1))},  축퇴 {2*l+1}")
```
![구면 위 입자의 축퇴](/assets/image-67.png)
```
l = 0: E =   0.0000 eV,  m = [0],  축퇴 1
l = 1: E =   7.6200 eV,  m = [-1, 0, 1],  축퇴 3
l = 2: E =  22.8599 eV,  m = [-2, -1, 0, 1, 2],  축퇴 5
l = 3: E =  45.7198 eV,  m = [-3, -2, -1, 0, 1, 2, 3],  축퇴 7
```
 
축퇴가 생기는 이유는 4.3절과 같습니다. 대칭성 때문이죠. 구는 어느 축으로 돌려도 같으니 $z$축을 어디로 잡든 에너지가 달라질 이유가 없습니다. 
 
그러니 대칭성을 깨면 축퇴도 깨집니다. 자기장을 걸면 특정 방향이 특별해지고, 그러면 $m$이 다른 상태들의 에너지가 갈라지죠. 이것을 **Zeeman 효과**라고 부릅니다. 그래서 $m$을 **자기양자수(magnetic quantum number)**라고 부르게 됩니다.


## 각운동량이 두 개
 
각운동량의 크기의 제곱에 해당하는 연산자 $\hat{L}^2$를 구면좌표로 쓰면 이렇게 됩니다.
 
$$
\hat{L}^2 = -\hbar^2\Lambda^2
$$
 
그러니 우리 방정식은 사실
 
$$
\hat{H} = \frac{\hat{L}^2}{2I}
$$
 
였던 셈이죠. 고전역학에서 회전 운동에너지가 $L^2/2I$인 것과 똑같습니다. 그리고 구면조화함수는 $\hat{L}^2$의 고유함수입니다.
 
$$
\hat{L}^2\,Y_l^m = l(l+1)\hbar^2\,Y_l^m
$$
 
각운동량의 크기가 $\sqrt{l(l+1)}\,\hbar$라는 뜻입니다. 한편 $\hat{L}_z$는 앞 절에서 본 그대로이고
 
$$
\hat{L}_z\,Y_l^m = m\hbar\,Y_l^m
$$
 
구면조화함수는 두 연산자의 공통 고유함수입니다. 실제로 $[\hat{L}^2, \hat{L}_z] = 0$라서 각운동량의 크기와 $z$ 성분을 동시에 알 수 있죠. 양자수가 $l$과 $m$ 두 개인 이유가 이것입니다.
 
그런데 세 축을 다 알 수는 없습니다. $\hat{L}_x$, $\hat{L}_y$, $\hat{L}_z$는 서로 교환하지 않거든요. 예를 들면, 이렇게요.
 
$$
[\hat{L}_x, \hat{L}_y] = i\hbar\hat{L}_z
$$
 
$z$ 성분을 확정하면 $x$와 $y$ 성분은 흐려집니다. 각운동량 벡터가 $z$축 주위 어딘가를 가리키고 있다는 것만 알 뿐, 정확히 어느 방향인지는 알 수 없죠.


## 다음 이야기

분자의 진동 운동이 적외선 분광학과 이어진 것처럼 회전 운동도 마이크로파 분광학과 이어져 있습니다. 그리고 기체 상태에서 적외선 분광 스펙트럼을 찍으면 진동과 회전이 함께 바뀌는 전이가 나타나 피크가 여러 갈래로 갈라지는 경우도 있습니다. 다음 절에서는 이에 대해 살펴보겠습니다.


## 확인 문제
 
1. $l=3$인 상태는 몇 개일까요? 각각의 $m$ 값을 나열해보세요.
2. 구면조화함수가 직교한다는 것을 수치적분으로 확인해보세요. (주의: 구면 위에서 적분하려면 $\sin\theta \, d\theta \, d\phi$로 적분해야 합니다.)
3. 구면조화함수의 마디는 어떻게 셀까요? $Y_2^0$과 $Y_2^2$를 그려보고 마디면의 개수와 방향을 비교해보세요.
