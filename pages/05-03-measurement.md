# 5.3. 측정에서 얻는 것

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/05-03.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞의 두 절에서 도구를 갖췄습니다. 연산자가 무엇인지, 왜 에르미트여야 하는지 알았죠. 이제 이 장의 진짜 질문에 답할 차례입니다.

> 실제로 측정하면 무슨 값이 나오는가?

2.3절에서 우리는 "전자가 어디 있는가"라는 질문을 포기했었죠. 하지만 측정을 하면 위치를 알 수 있습니다. 그러면 대체 입자의 위치나 운동량을 측정했을 때 우리가 얻는 값은 어디서 나오는 걸까요?


## 고유상태와 측정

가장 단순한 경우부터 봅시다. 계가 어떤 연산자($\hat{A}$)의 고유상태($\psi$)에 있는 경우입니다.

$$
\hat{A}\psi = a\psi
$$

이때 $\hat{A}$에 해당하는 물리량을 측정하면 반드시 $a$가 나옵니다. 몇 번을 측정해도 똑같이 $a$입니다.  

5.1절에서 상자 속 입자의 바닥상태에 $\hat{T}$를 적용했더니 원래 함수의 0.376 eV 배가 나왔던 것 기억하시나요? 5.2절에서 에르미트 연산자의 고유값이 실수라는 것을 확인했었죠. 상자 안에서는 $\hat{H} = \hat{T}$이니 이 상태의 에너지를 재면 항상 0.376 eV가 나옵니다. 이제 앞 절의 이야기가 왜 필요했는지 알겠습니다.

여기까지는 별로 이상한 건 없는 것 같네요. 문제는 다음입니다.


## 고유상태가 아니라면?

상자 속 입자의 바닥상태는 $\hat{H}$의 고유상태였습니다. 그런데 $\hat{x}$의 고유상태는 아니었죠. 5.1절에서 $x\psi$를 그렸을 때 모양이 바뀌었으니까요. 그러면 이 상태에서 위치를 측정하면 무엇이 나올까요?

사실 여기서 우리가 말할 수 있는 것은 어떤 값이 나올 확률이 얼마인가뿐입니다. 이 상황을 다루려면 5.2절의 마지막 결과가 필요합니다. 에르미트 연산자의 고유함수들은 서로 직교하고 완전기저를 이루니, 어떤 상태든 그것으로 펼칠 수 있습니다. 이렇게 고유상태들의 선형결합으로 이루어진 상태를 **중첩 상태(superposition state)**라고 부릅니다.

$$
\psi = \sum_n c_n\psi_n \qquad c_n = \langle\psi_n|\psi\rangle
$$

그리고 여기서 **양자역학의 핵심 가정** 중 하나가 등장합니다.

> $\hat{A}$에 해당하는 물리량을 측정하면 고유값 $a_n$ 중 하나가 나온다. 그리고 특정 $a_n$이 나올 확률은 $|c_n|^2$이다.

2.3절에서 $|\psi(x)|^2$이 위치의 확률밀도라고 했던 것을 기억하시나요? 그것은 위치 연산자에 대한 특수한 경우였고, 이것은 모든 물리량으로 확장된 형태입니다.

$$
\sum_n |c_n|^2 = 1
$$

그리고 당연하게도, 각각의 확률을 모두 더하면 1이 되어야 합니다. 


## 중첩 상태 만들기

상자 속 입자의 중첩 상태를 하나 만들어보죠. 바닥 상태와 첫 번째 들뜬 상태를 섞어보겠습니다.

$$
\psi = c_1\psi_1 + c_2\psi_2
$$

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e, h

L, N = 1e-9, 2000
x = np.linspace(0, L, N)


def psi(n):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)


fig, axes = plt.subplots(1, 3, figsize=(13, 4))

for ax, w in zip(axes, [1.0, 0.5, 0.25]):
    a, b = np.sqrt(w), np.sqrt(1 - w)
    P = a * psi(1) + b * psi(2)
    ax.plot(x * 1e9, P / np.abs(P).max(), color="crimson", lw=2, label="psi")
    ax.fill_between(
        x * 1e9, (P**2) / (P**2).max(), color="steelblue", alpha=0.35, label="|psi|^2"
    )
    ax.axhline(0, color="black", lw=0.6)
    ax.set_xlabel("x (nm)")
    ax.set_yticks([])
    ax.set_title(f"|c1|^2={w:.2f}, |c2|^2={1-w:.2f}", fontsize=10)

axes[0].legend(fontsize=8)
plt.show()
```
![중첩 상태](/assets/image-54.png)

이 중첩 상태들의 에너지를 측정하면 서로 다른 값들이 나올 겁니다. 왼쪽 상태는 고유상태 하나만 있으니 위에서 말한 것처럼 항상 0.376 eV로 나올 거고요. 하지만 다른 두 상태는 중첩 상태이니 각각의 에너지가 확률적으로 측정될 것입니다. 예를 들어, 맨 오른쪽 상태의 에너지를 측정하면 0.376 eV가 75%, 1.504 eV가 25%의 확률로 나온다는 뜻이죠. 그리고 0.376과 1.504 사이의 값은 절대 나오지 않습니다. 측정값은 언제나 고유값 중 하나만으로 나옵니다.


## 기댓값의 계산

그런데 우리가 실험에서 얻는 것은 보통 여러 번 측정한 평균입니다. 분광기로 스펙트럼을 찍으면 수많은 분자에서 나온 신호가 겹쳐 있죠. 이 값을 **기댓값(expectation value)**이라고 부르고 이렇게 씁니다.

$$
\langle A\rangle = \langle\psi|\hat{A}|\psi\rangle = \int\psi^*\hat{A}\psi\,dx
$$

중첩 상태로 계산해보면 왜 이것이 평균인지 드러납니다. $\psi = \sum_n c_n\psi_n$을 넣고 직교성을 쓰면

$$
\langle A\rangle = \sum_n |c_n|^2 a_n
$$

각 고유값에 그 확률을 곱해 더한 것입니다. 확률론에서 배운 평균의 정의 그대로죠.

시뮬레이션으로 확인해봅시다. 같은 상태를 20000번 측정한다고 하고, 각 측정에서 75%/25% 확률로 두 값 중 하나가 나오게 하면

```python
def En(n):
    return n**2 * h**2 / (8 * m_e * L**2) / e  # 4.1절의 결과, eV 단위


rng = np.random.default_rng(1)
draws = rng.choice([En(1), En(2)], size=20000, p=[0.75, 0.25])

fig, ax = plt.subplots(figsize=(7, 4))
ax.hist(draws, bins=40, color="steelblue", alpha=0.75)
ax.axvline(
    draws.mean(),
    color="crimson",
    lw=2.5,
    label=f"measured mean = {draws.mean():.3f} eV",
)
ax.axvline(
    0.75 * En(1) + 0.25 * En(2),
    color="black",
    lw=1.5,
    ls="--",
    label=f"<E> = {0.75 * En(1) + 0.25 * En(2):.3f} eV",
)
ax.set_xlabel("measured energy (eV)")
ax.set_ylabel("counts")
ax.legend(fontsize=9)
plt.show()

print(f"측정 평균 = {draws.mean():.4f} eV")
print(f"기댓값    = {0.75 * En(1) + 0.25 * En(2):.4f} eV")
```
![에너지 측정값의 평균과 기댓값](/assets/image-55.png)
```
측정 평균 = 0.6598 eV
기댓값    = 0.6581 eV
```

막대는 두 자리에만 서 있습니다. 그 중간에 있는 평균값이나 기댓값은 한 번도 측정되지 않는 값입니다. 주사위를 던져 나오는 눈의 기댓값은 3.5이지만 실제로 주사위를 던져서 3.5가 나오지는 않는 것과 같습니다. 기댓값은 "가장 나올 법한 값"이 아니라 "여러 번 재서 평균 낸 값"입니다.

측정을 20000번만 했으니 평균값이 기댓값과 조금 다르게 나오지만 측정 횟수를 더 크게 할수록 점점 더 가까워질 것입니다. 확률론에서 배웠던 것처럼요.


## 상자 속 입자의 기댓값들

이제 실제로 계산해봅시다. 고유상태에서 위치와 운동량의 기댓값이 얼마인지 봅시다.

```python
dx = L / (N + 1)
xg = np.linspace(0, L, N + 2)[1:-1]
D1 = (np.eye(N, k=1) - np.eye(N, k=-1)) / (2 * dx)
D2 = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2


def psi_g(n):
    return np.sqrt(2 / L) * np.sin(n * np.pi * xg / L)


unit = hbar / L  # 운동량의 자연 단위

for n in [1, 2, 3]:
    p = psi_g(n)
    x_mean = np.trapezoid(p * xg * p, xg)
    x2_mean = np.trapezoid(p * xg**2 * p, xg)
    p_mean = np.trapezoid(p * (-1j * hbar * (D1 @ p)), xg)
    p2_mean = np.trapezoid(p * (-(hbar**2) * (D2 @ p)), xg)
    print(
        f"n={n}:  <x>={x_mean/L:.4f} L   <x^2>={x2_mean/L**2:.4f} L^2   "
        f"<p>={abs(p_mean)/unit:.2e} (h/L)   <p^2>={p2_mean/unit**2:8.4f} (h/L)^2"
    )
```
```
n=1:  <x>=0.5000 L   <x^2>=0.2827 L^2   <p>=0.00e+00 (h/L)   <p^2>=  9.8696 (h/L)^2
n=2:  <x>=0.5000 L   <x^2>=0.3207 L^2   <p>=4.35e-16 (h/L)   <p^2>= 39.4784 (h/L)^2
n=3:  <x>=0.5000 L   <x^2>=0.3277 L^2   <p>=2.18e-16 (h/L)   <p^2>= 88.8263 (h/L)^2
```

두 결과가 눈에 띕니다.

$\langle x\rangle = L/2$입니다. 생각해보면 당연합니다. $|\psi_n|^2$이 가운데를 기준으로 좌우대칭이니까요. 하지만 실제로 전자가 가운데에 있다는 뜻은 아닙니다. $n=2$를 보세요. 가운데는 마디라 전자를 절대 발견할 수 없는 지점인데도 평균은 가운데입니다. 기댓값이 실제 측정값과 다를 수 있다는 좋은 예죠.

$\langle p\rangle = 0$입니다. 마찬가지로 대칭성 때문입니다. 전자가 오른쪽으로 갈 확률과 왼쪽으로 갈 확률이 같으니 평균 운동량이 0이죠. 상자 안에 갇혀 있으니 어디로도 흘러가지 않는다는 뜻입니다.

그런데 $\langle p^2\rangle$은 0이 아닙니다. 같은 단위로 쟀는데 $9.87$이 나왔죠. $10^{-16}$과 비교하면 하늘과 땅 차이입니다. 제곱하면 방향이 사라지니 왼쪽으로 가든 오른쪽으로 가든 똑같이 더해지기 때문입니다. 숫자를 자세히 보면 $9.8696 = \pi^2$이고, $39.4784 = 4\pi^2$, $88.8263 = 9\pi^2$입니다. 즉
 
$$
\langle p^2\rangle = \left(\frac{n\pi\hbar}{L}\right)^2
$$
 
de Broglie 관계 $p = \hbar k$에 $k_n = n\pi/L$을 넣은 것과 정확히 같죠. 그리고 이것을 $2m$으로 나누면
 
$$
\langle T\rangle = \frac{\langle p^2\rangle}{2m} = \frac{n^2h^2}{8mL^2}
$$
 
4.1절에서 손으로 구한 에너지 그대로입니다.

마지막으로 중첩 상태의 에너지 기댓값을 수치적분으로 직접 구해보죠. 앞에서 가중평균으로 얻은 값과 맞는지 확인하는 겁니다.
 
```python
c1, c2 = np.sqrt(0.75), np.sqrt(0.25)
Psi_g = c1 * psi_g(1) + c2 * psi_g(2)
H = -(hbar**2) / (2 * m_e) * D2

E_mean = np.trapezoid(Psi_g * (H @ Psi_g), xg) / e

print(f"수치적분으로 구한 <E>           = {E_mean:.5f} eV")
print(f"가중평균 0.75 * E1 + 0.25 * E2 = {0.75 * En(1) + 0.25 * En(2):.5f} eV")
```
```
수치적분으로 구한 <E>           = 0.65805 eV
가중평균 0.75 * E1 + 0.25 * E2 = 0.65805 eV
```
 
정확히 일치합니다. $\langle A\rangle = \sum_n|c_n|^2a_n$이라는 식이 실제로 성립하는 것이죠.
 
그런데 이 결과를 보면서 한 가지 눈여겨볼 것이 있습니다. $\langle E\rangle = 0.658$ eV는 바닥상태 에너지 0.376 eV보다 큽니다. $E_1$이 가장 작은 고유값이니 어떻게 섞든 평균이 그보다 작아질 수는 없거든요.
 
$$
\langle E\rangle = \sum_n |c_n|^2E_n \ge \sum_n |c_n|^2E_1 = E_1
$$
 
아무 함수나 가져와 에너지 기댓값을 계산해도 이 부등식이 성립합니다. 어떤 함수든 고유함수들로 펼칠 수 있으니까요. 이 성질을 **변분원리(variational principle)**라고 부르는데, 정답을 몰라도 여러 시험 함수를 넣어보고 가장 낮은 에너지를 주는 것을 고르면 된다는 뜻이라 앞으로 나올 근사 계산의 엔진이 됩니다. 11장에서 제대로 다루겠습니다.


## 다음 이야기

$\langle p\rangle = 0$인데 $\langle p^2\rangle \neq 0$이라는 것이 무슨 뜻일까요? 통계 과목에서 배웠던 **분산(variance)**이라는 개념을 생각해보면, 분산이 0이 아니라는 얘기는 편차가 존재한다는 뜻입니다. 분산을 어떻게 계산했었는지 기억하시나요? 제곱한 값의 평균과 평균의 제곱의 차이였습니다. 식으로 나타내면 이렇게 되죠: $\mathbf{Var}(\hat{p}) = \langle p^2\rangle - \langle p\rangle^2$. 통계에서는 단순히 데이터가 얼마나 퍼져있냐를 따지는 개념이지만 양자역학에서는 더 심오한 결과를 낳게 됩니다.


## 확인 문제

1. $\psi = \sqrt{0.5}\,\psi_1 + \sqrt{0.5}\,\psi_3$인 상태의 에너지 기댓값을 계산해보세요. 어떤 값들이 측정되고 확률은 각각 얼마인가요?
2. $n=2$ 상태에서 전자를 상자 왼쪽 절반에서 발견할 확률을 계산해보세요. $\langle x\rangle = L/2$라는 사실과 어떻게 어울리나요?
3. 두 고유상태를 섞은 상태의 에너지 기댓값이 결코 $E_1$보다 작아질 수 없는 이유를 $\sum|c_n|^2 = 1$을 이용해 설명해보세요.
4. (선택) 사인 함수가 아닌 아무 함수, 예를 들어 $f = x(L-x)$를 규격화해서 에너지 기댓값을 계산해보세요. $E_1$보다 큰가요?