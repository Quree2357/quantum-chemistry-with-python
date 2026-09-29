# 14.1. 오비탈이 섞이는 조건

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/14-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

12장에서 분자 오비탈을 구할 때 적분이 세 개가 등장했었습니다. Coulomb 적분 $H_{AA}$, 공명 적분 $H_{AB}$, 그리고 중첩 적분 $S$였죠. 그때는 이 적분값들을 정확하게 계산했었지만 분자가 커지면 계산량이 너무 늘어납니다. 그래서 이번 절에서는 적분값을 직접 계산하는 대신 매개변수로 두고, 이것들이 어떤 역할을 하는지부터 살펴보겠습니다.


## 적분 삼형제, 다시 돌아오다

세 개의 적분의 정의를 다시 살펴봅시다. 그리고 새로운 기호를 하나씩 붙여주도록 하죠.

$$
\alpha_A = H_{AA} = \left\langle \phi_A \middle| \hat{H} \middle| \phi_A \right\rangle, \quad \alpha_B = H_{BB} = \left\langle \phi_B \middle| \hat{H} \middle| \phi_B \right\rangle
$$
$$
\beta = H_{AB} = \left\langle \phi_A \middle| \hat{H} \middle| \phi_B \right\rangle, \quad S = \left\langle \phi_A \middle| \phi_B \right\rangle
$$

먼저 $\alpha$를 봅시다. 이 적분은 잘 보면 Hamiltonian의 기댓값을 계산하는 것과 같은 형태입니다. 즉, 분자 안에서 원자 오비탈에 있는 전자가 가지는 대략적인 에너지라고 볼 수 있습니다. 오비탈의 에너지라고도 볼 수 있죠. 실제로 계산해보면 그 의미가 더욱 분명하게 드러납니다. 12.2절에서 봤던 수소 분자 이온에 대해 다시 계산해보죠. 여기서는 $H_{AA}=H_{BB}$였습니다.
```python
import numpy as np


def S(R):
    return np.exp(-R) * (1 + R + R**2 / 3)


def H_AA(R):
    return -0.5 + np.exp(-2 * R) * (1 + 1 / R)


def H_AB(R):
    return -0.5 * S(R) - np.exp(-R) * (1 + R) + S(R) / R


print(f"{'R':>6} {'alpha':>8} {'beta':>8} {'S':>8}")
for R in [1.0, 2.0, 2.4935, 3.0, 4.0, 8.0]:
    print(f"{R:6.2f} {H_AA(R):8.4f} {H_AB(R):8.4f} {S(R):8.4f}")
```
```
     R    alpha     beta        S
  1.00  -0.2293  -0.3066   0.8584
  2.00  -0.4725  -0.4060   0.5865
  2.49  -0.4904  -0.3341   0.4599
  3.00  -0.4967  -0.2572   0.3485
  4.00  -0.4996  -0.1389   0.1893
  8.00  -0.5000  -0.0068   0.0102
```
$R$이 커짐에 따라 $\alpha$의 값이 빠르게 -0.5로 수렴하고 있습니다. 에너지를 계산했던 과정을 생각해보면 이 값은 수소 원자의 $1s$ 오비탈의 에너지와 같음을 알 수 있죠.

$$
E = \frac{H_{AA} \pm H_{AB}}{1 \pm S} \quad E_{R\to\infty} \approx H_{AA}
$$

12.2절에서 봤던 것처럼 이 값은 전자가 무한대의 거리에 있을 때에 비해 얼마나 에너지가 감소했는지에 대한 값입니다. 그러니까 부호만 반대고 이온화 에너지와 같은 크기입니다. 사실은 $R$이 무한히 커져야 정확히 같은 크기를 가지겠지만 결합 거리($R=2.49$) 근처에서도 충분히 비슷하니 어느 정도 근사값으로 써도 됩니다. 수소 원자뿐만 아니라 다른 원자에서도 마찬가지로요.  

그리고 $\beta$와 $S$는 $R$이 커지면 0으로 수렴합니다. $\beta$는 두 원자 오비탈 사이의 상호작용에 대한 값이라고 볼 수 있죠. 일반적인 범위에서는 음수 값을 가집니다. 이 값이 0에 가까워진다는 것은 Hamiltonian을 통해 두 오비탈이 거의 섞이지 않는다는 뜻입니다. $S$는 예전에 설명했던 것처럼 두 오비탈이 공간적으로 얼마나 겹쳐있는지를 나타내는 값이고요. 그래서 이 두 값은 $R$이 너무 작지만 않으면 대체로 비슷한 증감을 나타내는 경향을 보입니다. 실제로도 계산해보면 두 값의 비($\beta/S$)가 약 -0.7 근처로 나타납니다.  

우리는 여기서 계산을 단순화하기 위해 $S=0$으로 둘 수 있습니다. *엥? 오비탈이 안 겹친다고요? 그게 돼요?* 물론 그냥 띡하고 0으로 둘 수 있는 건 아니고요. 원자 오비탈 간의 중첩 효과를 공명 적분 $\beta$ 안으로 포함시키면 그렇게 근사시킬 수 있습니다. 이렇게 놓는 이유는 영년 방정식을 간단하게 만들어서 훨씬 다루기 쉽게 하려고입니다. 직접 예시를 살펴보면서 확인해봅시다.


## 이세계 트럭 대신 이원자 분자에 치이기

수소 분자처럼 동핵 이원자 분자이면 $A=B$라서 $H_{AA}=H_{BB}$였는데 HF나 CO처럼 두 원자가 다른 원자면 두 값이 달라질 겁니다. 일반성을 잃지 않고 원자 $A$의 오비탈 에너지가 더 낮다고 합시다.

$$
\Delta \alpha = \alpha_B - \alpha_A > 0
$$

여기서 $S=0$으로 두면 영년 방정식은 이렇게 되죠.

$$
\det \left( \begin{bmatrix}
\alpha_A - E & \beta \\
\beta & \alpha_B - E
\end{bmatrix}
\right) = (\alpha_A - E)(\alpha_B - E) - \beta^2 = 0
$$

마찬가지로 $E$에 대해 풀면 두 개의 해가 나옵니다.

$$
E_{\pm} = \frac{\alpha_A + \alpha_B}{2} \mp \sqrt{\left( \frac{\Delta \alpha}{2} \right)^2 + \beta^2}
$$

우리가 지금 관심 있는 것은 결합성 오비탈의 에너지가 원래보다 얼마나 더 내려가느냐입니다. 이 안정화 에너지를 $\delta$라고 하면 다음과 같이 쓸 수 있습니다.

$$
\delta = \alpha_A - E_+ = \sqrt{\left( \frac{\Delta \alpha}{2} \right)^2 + \beta^2} - \frac{\Delta \alpha}{2}
$$

극단적인 두 경우를 살펴봅시다. 먼저 $\Delta \alpha$가 아주 작으면 $\delta \approx |\beta|$가 되죠. 동핵 이원자 분자처럼 되는 거죠. 그리고 반대로 $\Delta \alpha$가 아주 커지면 제곱근의 이항 전개로부터 $\delta \approx \frac{\beta^2}{\Delta\alpha}$가 됩니다. 이 모양, 어디선가 본 것 같지 않나요?  
11.3절에서 봤던 2차 섭동 에너지와 똑같은 모양입니다.

$$
E_n^{(2)} = \sum_{m \neq n} \frac{\left| \hat{H}_{mn} \right|^2}{E_n^{(0)} - E_m^{(0)}}
$$

두 원자 오비탈 사이의 에너지 차이가 클 때는 $\beta$를 섭동으로 취급할 수 있는 거죠. 즉, 오비탈 간의 에너지 차이가 작을수록 잘 섞인다는 뜻이 됩니다.

LCAO의 계수도 확인해볼까요? 두 계수 간의 비율은 영년 방정식의 아무 줄에다 에너지 값을 대입하면 됩니다.

$$
\phi = c_A \phi_A + c_B \phi_B
$$
$$
\frac{c_B}{c_A} = \frac{\delta}{|\beta|}
$$

$\Delta \alpha=0$이면 $\delta=|\beta|$이니 $c_A=c_B$가 되어 두 오비탈이 절반씩 섞입니다. 그리고 $\Delta \alpha$가 커질수록 $\delta$가 작아지니 $c_B$의 값이 점점 작아지고 결합성 오비탈에는 $A$가 더 많이 기여하게 되죠. 반결합성 오비탈은 반대로 $B$가 더 많이 기여하게 되고요.

직접 계산해봅시다. 에너지 단위를 $|\beta|$로 잡고 $\alpha_A=0,\,\alpha_B=\Delta\alpha$로 두면 됩니다.
```python
import matplotlib.pyplot as plt


def orbital_mix(a_A, a_B, beta):
    H = np.array([[a_A, beta], [beta, a_B]])
    E, C = np.linalg.eigh(H)
    return E, C


Delta_a = np.linspace(0, 10, 200)
beta = -1.0
delta, c_A_squared = [], []
for x in Delta_a:
    E, C = orbital_mix(0, x, beta)
    delta.append(-E[0])
    c_A_squared.append(C[0, 0] ** 2)


fig, axes = plt.subplots(1, 2, figsize=(11, 5))

ax = axes[0]
ax.plot(Delta_a, delta, lw=3, color="steelblue", label="exact")
ax.plot(Delta_a[10:], beta**2 / Delta_a[10:], ls="--", lw=2, color="darkorange", label="β^2 / Δα")
ax.set_ylim(0, 1.1)
ax.set_xlabel("Δα / |β|")
ax.set_ylabel("δ / |β|")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(Delta_a, c_A_squared, lw=3, color="steelblue")
ax.axhline(0.5, color="gray", lw=1)
ax.set_ylim(0.5, 1.02)
ax.set_xlabel("Δα / |β|")
ax.set_ylabel("weight of A (bonding)")
ax.grid(alpha=0.3)

plt.show()
```
![이핵 이원자 분자에서 오비탈의 섞임](/assets/image-98.png)

$\Delta \alpha$가 $|\beta|$의 두 배만 되어도 안정화 에너지가 $0.5|\beta|$ 이하로 줄어들고, 결합성 오비탈의 85% 이상이 원자 $A$가 기여합니다. 그리고 $\Delta \alpha$가 커질수록 섭동론의 근사식 $\frac{\beta^2}{\Delta\alpha}$과 비슷해지죠.

여기서 $\alpha$의 의미를 다시 한번 생각해봅시다. 이온화 에너지에 음수를 붙인 값이었죠. 그러니까 $\alpha$가 낮은 원자는 전자를 더 세게 잡고 있다는 뜻입니다. 이걸 일반화학 시간에 뭐라고 불렀는지 기억하시나요? 바로 전기음성도였습니다!  

그러니까 다시 말하면, $\alpha$가 낮은 원자가 결합성 오비탈에 기여하는 정도가 더 크다는 것이고, 결합에 참여하는 전자가 전기음성도가 큰 원자 쪽에서 더 많이 발견된다는 뜻입니다. 이게 바로 극성 공유 결합의 정체죠. 전기음성도 차이가 더욱 커지면 결합성 오비탈이 거의 전부 원자 하나의 오비탈이 되니 이온 결합이 되는 겁니다.

[[TIP]]
사실 첫 부분에도 말했지만 $\alpha$가 정확하게 이온화 에너지와 같은 값을 가지는 건 아닙니다. 그래서 전기음성도를 정의할 때는 이온화 에너지만 가지고 정의하지는 않는데요. Pauling이 정의한 전기음성도는 결합 에너지로부터, Mulliken의 정의는 이온화 에너지와 전자 친화도를 같이 사용하죠.
[[/TIP]]


## 오비탈이 섞이는 조건

드디어 이 절의 제목에 주목해봅니다. 두 오비탈이 얼마나 잘 섞이는지를 판단하려면 안정화 에너지 $\delta$를 보면 되죠.

$$
\delta = \sqrt{\left( \frac{\Delta \alpha}{2} \right)^2 + \beta^2} - \frac{\Delta \alpha}{2}
$$

이 식에서 안정화 에너지가 커지려면 $\Delta \alpha$가 작아지거나 $|\beta|$가 커져야 합니다. $\beta$는 두 가지 요인에 의해서 결정되는 값이었죠. 중첩 적분 $S$의 값과 비슷한 경향을 가지지만 두 오비탈의 대칭 종이 다르면 정확히 0이 되기도 합니다. 그러니까 오비탈이 잘 섞이려면 다음 세 가지 조건을 만족해야 하겠군요.

1. 에너지가 비슷해야 한다: $\Delta \alpha$가 작아야 한다.
2. 중첩이 커야 한다: $S$가 커져야 $|\beta|$가 커진다.
3. 대칭 종이 같아야 한다: 그래야 $\beta \neq 0$이 된다.

첫 번째와 두 번째 조건은 정도의 문제지만 마지막 조건은 0이냐 아니냐를 결정해버립니다. 에너지가 아무리 비슷하고 중첩이 아무리 커도 대칭 종이 다르면 전혀 섞이지 않는 거죠. 그래서 실제로 분자 오비탈을 구성할 때는 먼저 대칭 종으로 블록을 나누고, 그다음에 각 블록 안에서 나머지 두 조건으로 얼마나 섞이는지를 봅니다.


## 다음 이야기

이제 오비탈이 언제 섞이고 얼마나 섞이는지를 판단해 볼 수 있는 도구를 얻었습니다. 다음 절에서는 이걸 사용해 2주기 원소들의 동핵 이원자 분자들을 살펴보겠습니다. 그리고 액체 산소가 왜 자석에 끌리는지도 알 수 있죠.


## 확인 문제

1. 반결합성 오비탈에서 원자 $B$가 차지하는 비중을 계산해보세요. 결합성 오비탈에서 원자 $A$가 차지하는 비중과 비교하면 어떤가요?
2. $S \neq 0$일 때도 결합성 오비탈이 $\alpha$가 낮은 원자 쪽으로 쏠리는지 SciPy의 `eigh(H, S)`로 확인해보세요. $S=0.2$, $\beta=-1$로 두고 $\Delta \alpha$를 바꿔가며 계산하면 됩니다.
3. 수소의 이온화 에너지는 13.60 eV, 리튬은 5.39 eV, 플루오린은 17.42 eV입니다. $|\beta|$ 값이 비슷하다고 가정하면 결합성 오비탈이 한쪽으로 쏠리는 정도가 LiH와 HF 중 어느 분자에서 더 클까요? 각각 어느 원자 쪽으로 쏠릴까요?
