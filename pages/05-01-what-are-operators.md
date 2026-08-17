# 5.1. 연산자란 무엇인가

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/05-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

0.2절 도입부에서 연산자(operator)는 함수를 다른 함수로 바꿔주는 것이라고 했었습니다. 그리고 선형 연산자는 행렬로 쓸 수 있기 때문에 미분 연산자를 행렬로 써봤고요. Schrödinger 방정식에서는 Hamiltonian이라는 연산자가 등장했었습니다. 그런데 연산자는 사실 단순한 함수 변환기 그 이상의 의미를 가지고 있습니다. 이번에는 연산자에 대해서 더 자세히 알아보도록 하죠.


## 연산자와 선형성

원래 수학에서 연산자라고 하는 것은 반드시 함수에만 적용되는 것은 아닙니다. 선형대수학에서 들어보셨을 벡터 공간에서의 사상(map)이라는 것도 연산자의 일종입니다. 하지만 물리학에서는 연산자는 대부분 함수에 적용되는 경우가 많습니다. 특히 양자역학에서는 더욱 중요하고요. 연산자에는 연산자라는 표시를 위해 모자 기호(^)를 위에 붙인다고 했었죠.

$$
\hat{A}f(x) = g(x)
$$

이런 형태로 말이죠.  

양자역학에서는 선형 연산자(linear operator)만 다루게 됩니다. 연산자가 선형이라는 것의 정확한 정의는 연산자가 다음 조건을 만족할 때입니다.

$$
\hat{A} \left( af(x) + bg(x) \right) = a\hat{A}f(x) + b\hat{A}g(x)
$$

선형대수학에서의 선형 사상과 정확히 같은 정의입니다. 친숙하시죠?  


## 양자역학에서의 연산자

양자역학에서 가장 자주 쓰이는 연산자 세 가지는 위치, 운동량, 그리고 에너지 연산자입니다. 3차원 공간에서라면 원래는 각 축마다 하나씩 나오겠지만 모양은 똑같으니 $x$ 방향만 살펴보죠.

### 위치 연산자

별거 없습니다. 그냥 파동함수에 $x$를 곱해줍니다.

$$
\hat{x} = x \quad \Longrightarrow \quad \hat{x}\psi = x\psi
$$

너무 당연해 보여서 굳이 연산자라고 해야 하나 싶지만... 함수($\psi$)를 다른 함수($x\psi$)로 바꿔주니 연산자가 맞긴 합니다.

### 운동량 연산자

고전적인 운동량은 $p_x = m v_x$였었죠. 양자역학에서는 운동량에 대응되는 연산자가 모양이 조금 다릅니다.

$$
\hat{p} = -i\hbar\frac{d}{dx}
$$

*엥? 운동량에 왜 $i$가 붙어요? 운동량은 실수잖아요!*  
2.1절에서 배웠던 de Broglie 관계식 $p=\hbar k$를 떠올려봅시다. 그리고 파수 $k$를 가지는 파동은 $e^{ikx}$였죠. 이 함수를 미분해보면 다음과 같습니다.

$$
\frac{d}{dx}e^{ikx} = ike^{ikx}
$$

앞에 허수 단위 $i$가 딸려 나왔죠. 그래서 앞에 $-i$를 곱해줍니다. 그리고 $\hbar$를 곱해주면,

$$
-i\hbar\frac{d}{dx}e^{ikx} = \hbar ke^{ikx} = pe^{ikx}
$$

깔끔하게 앞에 $p$가 나옵니다. 즉, 파동에서 운동량을 뽑아내려면 운동량 연산자는 $-i\hbar$가 필요했던 겁니다.  

### 에너지 연산자

운동량을 알면 운동 에너지는 자동으로 나옵니다. 퍼텐셜 에너지는 우리가 $V(x)$로 주니까 이미 알고 있는 거고요.  
고전물리학에서 $T=\frac{1}{2}mv^2=\frac{p^2}{2m}$이었으니 이걸 똑같이 쓰면 됩니다.

$$
\hat{T} = \frac{1}{2m}\hat{p}\hat{p} = \frac{1}{2m}\left(-i\hbar\frac{d}{dx}\right)\left(-i\hbar\frac{d}{dx}\right) = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2}
$$

2차 미분이 나왔습니다. 그리고 여기에 퍼텐셜 에너지를 더해주면 Hamiltonian을 얻죠.

$$
\hat{H} = \hat{T} + \hat{V} = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)
$$

앗, 뭔가 규칙이 보이시지 않나요? 고전역학의 표현에서 $x \to \hat{x}$로, $p \to \hat{p}$로 바꾸면 양자역학의 연산자 표현이 되네요!  


## 직접 눈으로 확인해보기

위에서 본 연산자들이 구체적으로 어떤 역할을 하는지 직접 확인해봅시다. 4장에서 풀었던 상자 속 입자 모델에 한번 적용해보죠.

$$
\psi_1(x) = \sqrt{\frac{2}{L}}\sin\frac{\pi x}{L}
$$

그리고 0.2.3에서 만든 미분 행렬을 씁시다.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e

L, N = 1e-9, 2000
dx = L / (N + 1)
x = np.linspace(0, L, N + 2)[1:-1]


def D1(N, dx):
    return (np.eye(N, k=1) - np.eye(N, k=-1)) / (2 * dx)


def D2(N, dx):
    return (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2


psi = np.sqrt(2 / L) * np.sin(np.pi * x / L)
```
이제 세 연산자들을 각각 적용해봅시다.
```python
x_psi = x * psi
p_psi = -1j * hbar * (D1(N, dx) @ psi)
T_psi = -(hbar**2) / (2 * m_e) * (D2(N, dx) @ psi)

fig, axes = plt.subplots(1, 3, figsize=(13, 4))

axes[0].plot(x * 1e9, psi / psi.max(), color="black", lw=2, label="psi")
axes[0].plot(
    x * 1e9, x_psi / x_psi.max(), color="crimson", lw=2, label="x psi (scaled)"
    )
axes[0].set_title("position")

axes[1].plot(x * 1e9, psi / psi.max(), color="black", lw=2, label="psi")
axes[1].plot(
    x * 1e9, p_psi.imag / np.abs(p_psi).max(), color="steelblue", lw=2, label="Im(p psi) (scaled)"
)
axes[1].set_title("momentum")

axes[2].plot(x * 1e9, psi / psi.max(), color="black", lw=2, label="psi")
axes[2].plot(
    x * 1e9, T_psi / T_psi.max(), color="seagreen", lw=2, ls="--", label="T psi (scaled)"
)
axes[2].set_title("kinetic energy")

for ax in axes:
    ax.axhline(0, color="gray", lw=0.6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("x (nm)")
    ax.set_yticks([])
    ax.legend(fontsize=8)

plt.show()
```
![위치, 운동량, 에너지 연산자](/assets/image-52.png)

위치 연산자는 오른쪽으로 무게중심이 옮겨갔습니다. $x$를 곱했으니 $x$가 큰 쪽이 더 커진 거죠.  

운동량 연산자는 사인 함수가 코사인 함수가 되어 모양이 완전히 바뀌었습니다. 그리고 실수부가 0인 순허수 값으로 나타났습니다.  

운동에너지 연산자는 원래 함수와 정확히 똑같은 모양이 되었습니다! 물론 크기는 달라졌지만요. 그림에서는 일부러 모양이 같다는 것을 보기 위해 비율을 맞춰서 그렸습니다.  
이게 무슨 뜻이었죠? 바로 고유함수라는 뜻이었습니다. 행렬을 곱해도 방향이 안 바뀌는 고유벡터와 마찬가지로 연산자를 적용해도 모양이 바뀌지 않는 함수죠. 그런데 크기가 얼마나 달라진 걸까요?
```python
ratio = (T_psi / psi) / e
print(f"T psi / psi = {ratio.mean():.6f} eV (표준편차 {ratio.std():.2e})")

E1 = np.pi**2 * hbar**2 / (2 * m_e * L**2) / e
print(f"해석해 E1    = {E1:.6f} eV")
```
```
T psi / psi = 0.376030 eV (표준편차 2.59e-09)
해석해 E1    = 0.376030 eV
```
그리고 그 비율은 4.1절에서 손으로 구했던 바닥상태의 에너지와 정확히 같습니다. Schrödinger 방정식($\hat{H}\psi=E\psi$)을 고유값 문제로 보면 에너지가 고유값으로 나타난다고 했었습니다. 지금 그 고유값이 나온 겁니다. 상자 안에서는 $V=0$이어서 $H=T$가 되거든요.  

하지만 위치와 운동량은 고유함수가 아니었습니다. 이 차이가 무슨 일을 만드는지는 다음 절에서 보게 될 것입니다.


## 어어 새치기는 안 된다

연산자는 일반적인 곱셈과 다르게 순서를 바꾸면 다른 결과값을 줍니다. 행렬이랑 비슷하죠? 사실 앞에서 연산자를 행렬로 써봤으니 이미 예상하셨을 수도 있겠지만요...  
양자역학에서는 연산자의 순서를 바꿔서 곱셈을 하는 경우가 생각보다 많이 있어서 아예 이름을 붙여버렸습니다.

$$
[\hat{A}, \hat{B}] = \hat{A}\hat{B} - \hat{B}\hat{A}
$$

이것을 **교환자(commutator)**라고 부릅니다. 그러니까, 어떤 함수에 B 연산자를 먼저 적용하고 A 연산자를 적용한 결과에서 A 연산자를 먼저 적용하고 B 연산자를 적용한 결과를 뺀 거죠. 이 값이 모든 함수에 대해 0이라면 두 연산자는 **교환된다(commute)**고 말합니다.  


## 다음 이야기

양자역학에서의 연산자에 대해 알아보았습니다. 그런데 아무 연산자나 가져와도 선형이기만 하면 될까요? 연산자에 대체 무슨 의미가 있는 걸까요? 운동량 연산자는 순허수 값을 갖는데 이게 무슨 뜻일까요? 다음 절에서 이 질문들에 답해보겠습니다.


## 확인 문제
1. $\hat{p}$를 상자 속 입자의 $n=2$ 상태에 적용해보세요. 결과가 어떻게 나타나나요?
2. $\hat{T}$를 상자 속 입자의 $n=2, 3$ 상태에 적용하고 원래 함수와의 비율을 구해보세요. 4.1절에서 손으로 구했던 에너지 값과 같나요? (주의: 마디 근처에서는 값이 너무 작아 나눗셈이 불안정해질 수 있습니다. 충분히 큰 점들만 골라서 평균을 내보세요.)
3. 자유 전자의 파동을 나타내는 함수 $e^{ikx}$는 $\hat{p}$의 고유함수입니다. 직접 확인해보고 고유값이 얼마인지 계산해보세요.
4. 두 연산자 $\hat{A}=\frac{d}{dx}$와 $\hat{B}=x^2$에 대해 교환자를 구해보세요. 이 연산자들은 교환되나요?
