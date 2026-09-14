# 12.2. 수소 분자 이온

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/12-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

수소 분자 이온은 두 개의 양성자와 한 개의 전자로 이루어져 있는 가장 간단한 분자입니다. 약 105.7 pm의 결합 길이와 -0.6026 hartree의 바닥 상태 에너지를 가지고 있는 안정한 분자죠.

```python (나중에 코드 삭제해야 됨)
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(5, 4))
ax.plot(-1, 0, "o", ms=18, color="crimson")
ax.text(-1, -0.25, 'nucleus A (+e)', fontsize=10, ha='center')
ax.plot(1, 0, "o", ms=18, color="crimson")
ax.text(1, -0.25, 'nucleus B (+e)', fontsize=10, ha='center')
ax.plot(0.7, 1, "o", ms=12, color="steelblue")
ax.text(0.7, 1.2, 'electron', fontsize=10, ha='center')

ax.annotate("", xy=(0.7, 1), xytext=(-1, 0), arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))
ax.annotate("", xy=(0.7, 1), xytext=(1, 0), arrowprops=dict(arrowstyle='->', lw=2.5, color='black'))
ax.text(-0.25, 0.6, '$r_A$', fontsize=14)
ax.text(0.9, 0.55, '$r_B$', fontsize=14)
ax.annotate("", xy=(-1, 0), xytext=(1, 0), arrowprops=dict(arrowstyle='<->', lw=2, color='darkorange'))
ax.text(0, -0.2, "R", fontsize=14, color='darkorange')

ax.set_xlim(-1.5, 1.5)
ax.set_ylim(-0.5, 1.2)
ax.set_aspect('equal')
ax.axis('off')
plt.show()
```
![수소 분자 이온](/assets/image-90.png)

Born-Oppenheimer 근사를 사용하면 Hamiltonian을 이렇게 쓸 수 있습니다.

$$
\hat{H} = -\frac{1}{2} \nabla^2 - \frac{1}{r_A} - \frac{1}{r_B} + \frac{1}{R}
$$

이제 전자에 대한 Schrödinger 방정식을 풀면 됩니다.


## 정확히 풀리는 마지막 순간

놀랍게도 수소 분자 이온에 대한 Schrödinger 방정식은 정확하게 풀 수 있습니다. 삼체 문제이긴 하지만 원자핵의 위치가 고정되어 있으니 좌표계를 적당히 설정하면 변수분리가 되거든요. 바로 타원 좌표계(elliptic coordinate system)를 쓰는 겁니다. 정확히는 장축 타원 좌표계(prolate spheroidal coordinate system)를 쓰는데요.

$$
\mu = \frac{r_A + r_B}{R} \qquad \nu = \frac{r_A - r_B}{R}
$$

이렇게 2차원 타원 좌표계를 잡은 다음에 $z$축으로 돌려서 만드는 좌표계입니다. 돌아간 각도를 $\phi$라고 하고요. 두 원자핵을 타원의 초점으로 생각하고 계산하는 방법입니다. 이 방식으로 $R$에 따라 분자의 에너지가 어떻게 변하는지 그려보면 아래처럼 나옵니다.

![수소 분자 이온의 PES](https://thumb.wikimedia.org/wikipedia/commons/thumb/0/04/H2plus_figure_2.png/960px-H2plus_figure_2.png)  
[출처: https://en.wikipedia.org/wiki/Dihydrogen_cation]

방정식을 풀 때 자기 양자수 $m$은 에너지에서는 제곱 항으로만 등장합니다. 즉, $m=0$인 상태를 제외하고 다른 상태들은 2중으로 축퇴되어 있다는 것이죠. 7.1절에서 봤던 고리 위 입자에서와 비슷하게요.  
그래서 이 $m$의 절댓값에 따라 이름을 붙이는데요, $m=0$이면 $\sigma$, $m=\pm1$이면 $\pi$, $m=\pm2$이면 $\delta$ 등과 같이 씁니다. 그리고 또 파동 함수의 대칭성에 따라서도 분류해서 원점을 중심으로 대칭인 파동 함수는 아래첨자 $g$를 쓰고, 반대칭인 경우에는 아래첨자 $u$를 씁니다.  
이 표기법에 따르면 위 그래프에서 가장 아래에 있는 빨간색 선은 $1\sigma_g$ 상태에 해당하고 그 바로 위에 있는 녹색 점선은 $1\sigma_u$ 상태에 해당합니다. 뭔가 익숙한 듯한 느낌이 드시죠?

그런데 이 방법은 원자핵이 두 개인 분자에서만 통합니다. 당장 물 분자만 해도 원자가 3개이니 이 좌표계를 쓸 수 없죠. 확장성이 없다는 겁니다. 그래서 우리는 근사법을 도입해서 이 문제를 풀어볼 겁니다. 그러면 나중에 좀 더 큰 분자에 대해서도 같은 방법을 쓸 수 있으니까요.


## 분자 오비탈의 탄생

원자의 Schrödinger 방정식을 풀면 오비탈이 나왔으니 분자에 대해서도 똑같은 접근을 할 수 있지 않을까요? 그렇게 탄생한 것이 **분자 오비탈 이론(molecular orbital theory)**입니다. 이 이론은 1930년대 초반에 독일의 물리학자 Hund와 미국의 화학자 Mulliken에 의해 만들어졌습니다. 분자 오비탈 이론에서는 전자가 화학 결합 안에 존재하는 것이 아니라 전체 분자에 퍼져서 존재한다고 봅니다.  

그렇다면 분자 오비탈은 어떻게 생겼을까요? 어떤 수학적인 형태로 표현할 수 있을까요? 여기서는 11.2절에서 봤던 선형 변분법을 사용해보겠습니다. 시험 함수를 정해야 할텐데, 기저 함수를 무엇으로 정하면 좋을까요?  

극단적인 상황을 가정해봅시다. $r_A \ll r_B$라면 전자가 느끼는 상황은 그냥 $A$ 위치에 양성자가 하나 있는 수소 원자의 경우와 비슷할 것입니다. 반대 경우도 마찬가지고요. 그러니까 이 두 경우를 적절히 조합하면 되지 않을까요?

$$
\phi = c_A 1s_A+c_B 1s_B
$$

$1s_A$와 $1s_B$는 각각 $A$와 $B$를 중심으로 하는 수소 원자의 오비탈입니다. 이런 식으로 선형 변분 함수, 즉 분자 오비탈을 만드는 방법을 **원자 오비탈의 선형 결합(linear combination of atomic orbitals)**이라고 부르고 LCAO라고 줄여 씁니다.


## 변분법으로 풀어보기

선형 변분법을 푸는 방법은 바로 영년 방정식을 푸는 것이었습니다.

$$
\det \left( \mathbf{H} - E \mathbf{S} \right) = 0
$$
$$
\det \left( \begin{bmatrix}
H_{AA} & H_{AB} \\
H_{BA} & H_{BB}
\end{bmatrix} - E \begin{bmatrix}
1 & S \\
S & 1
\end{bmatrix}
\right) = 0
$$

$1s_A$와 $1s_B$는 규격화되어 있다고 가정하면 $\mathbf{S}$의 대각 성분이 1이 되죠. 하지만 둘이 직교한다는 보장은 없기 때문에 나머지 요소들은 $S$라는 값으로 두었습니다.  
여기서 주목할 점은 수소 분자 이온이 동핵 이원자 분자라는 것입니다. 이 말은 원자핵 $A$와 $B$를 바꿔도 달라지는 게 없는 대칭성을 가진다는 것이죠. 따라서 $H_{AA}=H_{BB}$이고 $H_{AB}=H_{BA}$입니다.

$$
\det \left( \begin{bmatrix}
H_{AA} - E & H_{AB} - ES \\
H_{AB} - ES & H_{AA} - E
\end{bmatrix}
\right) = 0
$$
$$
(H_{AA} - E)^2 - (H_{AB} - ES)^2 = 0
$$

위 식을 정리해서 $E$에 대해 풀면 다음과 같이 두 개의 해를 얻습니다.

$$
E_+ = \frac{H_{AA} + H_{AB}}{1 + S}, \quad E_- = \frac{H_{AA} - H_{AB}}{1 - S}
$$

계수 $c$를 구하기 위해 원래 식에 두 개의 해를 대입하면 $E_+$의 경우에는 $c_A = c_B$가, $E_-$의 경우에는 $c_A = -c_B$가 나옵니다. 규격화까지 마치고 나면 두 개의 분자 오비탈을 얻습니다.

$$
\phi_+ = \frac{1s_A + 1s_B}{\sqrt{2(1+S)}}, \quad \phi_- = \frac{1s_A - 1s_B}{\sqrt{2(1-S)}}
$$

이렇게 얻은 두 개의 분자 오비탈 $\phi_+$와 $\phi_-$는 각각 결합성 오비탈과 반결합성 오비탈이라고 부릅니다. 왜 이런 이름을 가지게 되었는지는 다음 절에서 살펴보도록 하겠습니다.


## 적분 삼형제

계산 과정에서 나왔던 세 개의 적분 $H_{AA}$, $H_{AB}$, $S$를 봅시다. 이 값들을 타원 좌표계를 쓰면 계산할 수 있는데, 중간 과정은 생략하고 결과만 보도록 하죠. 우선 $S$는 이렇게 나타납니다.

$$
S = e^{-R} \left( 1 + R + \frac{R^2}{3} \right)
$$

11.2절에서 $S$를 중첩 적분이라고 불렀던 것을 기억하시나요? $S$의 값은 $R=0$일 때는 1이고 $R \to \infty$이면 0으로 수렴합니다. 원자핵 사이의 거리가 0이면 두 원자 오비탈이 완전히 겹치고, 무한대로 가면 완전히 따로 떨어진다는 의미죠. 그래서 중첩 적분이라고 부른 것입니다.

나머지 두 적분은 각각 Coulomb 적분과 공명(교환) 적분이라고 부릅니다.

$$
H_{AA} = -\frac{1}{2} + e^{-2R} \left(1 + \frac{1}{R}\right), \quad H_{AB} = -\frac{S}{2} - e^{-R}(1 + R) + \frac{S}{R}
$$

이제 위에서 구했던 에너지 식에 이 적분값들을 대입합시다. 그리고 그림으로 그려보죠.

```python
import numpy as np
import matplotlib.pyplot as plt


def S(R):
    return np.exp(-R) * (1 + R + R**2 / 3)


def H_AA(R):
    return -0.5 + np.exp(-2 * R) * (1 + 1 / R)


def H_AB(R):
    return -0.5 * S(R) - np.exp(-R) * (1 + R) + S(R) / R


def E_plus(R):
    return (H_AA(R) + H_AB(R)) / (1 + S(R))


def E_minus(R):
    return (H_AA(R) - H_AB(R)) / (1 - S(R))


R = np.linspace(0.5, 10, 5000)
Eg, Eu = E_plus(R), E_minus(R)

fig, ax = plt.subplots(figsize=(7, 6))
ax.plot(R, Eg, lw=3, color="steelblue", label="bonding (σ_g)")
ax.plot(R, Eu, lw=3, color="crimson", label="antibonding (σ_u*)")
ax.axhline(-0.5, color="black", ls="--", lw=2, label="H + H+ (-0.5)")

i = np.argmin(Eg)
ax.plot(R[i], Eg[i], "o", ms=9, color="darkorange")
ax.set_xlim(0.5, 10)
ax.set_ylim(-0.6, 0.2)
ax.set_xlabel("R (Bohr)")
ax.set_ylabel("E (Hartree)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()


print(f"R_e = {R[i]:.4f} Bohr = {R[i] * 0.529177:.4f} Angstrom")
print(f"E = {Eg[i]:.6f} Hartree = {Eg[i] * 27.2114:.4f} eV")
print(f"D_e = {(-0.5 - Eg[i]) * 27.2114:.4f} eV")
```
![수소 분자 이온의 1σ 분자 오비탈](/assets/image-91.png)
```
R_e = 2.4935 Bohr = 1.3195 Angstrom
E = -0.564831 Hartree = -15.3698 eV
D_e = 1.7641 eV
```

파란색으로 나타낸 결합성 오비탈에 최소점이 있네요. 이 상태의 에너지가 분리 상태(수소 원자 + 수소 원자 이온)보다 더 낮고요. 반면 빨간색으로 나타낸 반결합성 오비탈에서는 항상 분리 상태보다 에너지가 높습니다. 그래서 여기에 전자가 있으면 분자가 오히려 불안정해지죠.


## 변분법, 제 점수는요

결합 길이와 에너지를 실제 값과 비교해봅시다. 처음에 수소 분자 이온의 결합 길이는 약 105.7 pm이고 바닥 상태 에너지는 -0.6026 hartree라고 했었죠. 계산으로 구한 값은 131.95 pm와 -0.5648 hartree입니다. 결합 길이는 약 25% 정도 길게 나왔고, 에너지는 약 37% 정도 작게 나왔네요.

결합이 존재하고, 결합성 오비탈과 반결합성 오비탈이 있다는 것, 그리고 대략적인 크기 정도는 맞지만 정량적으로 쓰기에는 오차가 너무 큽니다. 이렇게 어긋난 이유는 바로 처음에 시험 함수를 정할 때 기저 함수를 두 개만 썼기 때문입니다. 지금 우리는 '수소 원자의 오비탈' 두 개를 선형 결합해서 시험 함수를 만들었습니다. 실제 수소 분자 이온에서 정확히 맞을 리가 없죠.

매개변수를 추가하거나 다른 기저 함수를 추가하는 방법 등 이를 개선하는 방법은 많습니다. 이 이야기는 나중에 다시 만나보도록 합시다.


## 다음 이야기

다음 절에서는 이 근사 결과를 이용해 수소 분자 이온에서 결합이 왜 생기는지를 더 자세히 살펴보도록 하겠습니다. 그리고 교과서에서 설명하는 결합이 생기는 이유가 정확하지 않다는 것도 확인하게 될 겁니다.


## 확인 문제

1. $R=0$일 때와 $R \to \infty$일 때 Coulomb 적분 $H_{AA}$와 공명 적분 $H_{AB}$가 어떻게 되는지 확인해보세요. 물리적으로는 무슨 의미일까요?
2. $E_+$와 $E_-$ 중 어느 쪽이 분리 상태의 에너지에서 더 많이 벗어나 있나요? 반결합성이 결합성보다 더 많이 벗어난 이유를 설명해보세요.
3. 반결합성 오비탈에도 최소점이 정말 존재하지 않나요? 코드로 확인해보세요.
4. 헬륨 분자 이온에도 같은 방법을 쓸 수 있을까요? 무엇이 달라져야 할까요?
