# 5.2. 에르미트 연산자와 브라켓 표기법

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/05-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

양자역학에서의 연산자는 모두 선형이라고 했습니다. 그러면 선형인 연산자라면 모두 양자역학에서 의미가 있을까요?  
물리량을 나타내는 연산자에는 조건이 하나 있습니다. 생각해보면 아주 당연한 것입니다. 바로 우리가 **실험으로 측정하는 값은 실수로 나와야 한다**는 것이죠. 아니, 실험을 실수하라는 게 아니라 복소수처럼 허수 단위가 붙으면 안 된다는 뜻입니다. 그럼 연산자는 선형이어야 한다는 것 말고 또 무슨 조건을 가져야 할까요?


## 에르미트 연산자

0.1.4에서 에르미트 행렬을 배웠습니다. 켤레전치가 자기 자신($A^\dagger = A$)인 행렬이었죠. 그리고 중요한 성질 두 가지가 있었습니다. 고유값이 반드시 실수이고, 고유벡터는 서로 직교한다는 것이었습니다.  
연산자에도 같은 개념이 있습니다. 다음 관계가 모든 함수 $f$와 $g$에 대해 성립하면 $\hat{A}$를 **에르미트 연산자(Hermitian operator)**라고 부릅니다.

$$
\int f^*\hat{A} g \, dx = \int \left(\hat{A} f \right)^* g \, dx
$$

그리고 행렬에서와 똑같은 결론이 따라옵니다.

> 에르미트 연산자의 고유값은 반드시 실수이고, 서로 다른 고유값에 속한 고유함수는 직교한다.

그래서 양자역학에서 물리량에 대응하는 연산자는 선형이어야 할 뿐만 아니라, 에르미트 연산자로 표현되어야 합니다.


## 브라켓 표기법

양자역학에서는 위와 비슷한 적분식이 하도 많이 나와서 매번 귀찮은 적분과 $dx$ 기호를 쓰는 대신 아예 다음과 같이 줄여서 표기합니다.

$$
\langle f | g \rangle \equiv \int f^* g \, dx
$$

이 표기법은 Dirac이 만들어서 **Dirac 표기법** 또는 **브라켓 표기법(bra-ket notation)**이라고 부릅니다. 이름의 유래는 예상하셨겠지만, 괄호를 뜻하는 영어 단어 bracket입니다. 이 단어를 반으로 쪼개서, 왼쪽 $\langle f|$를 **브라(bra)**, 오른쪽 $|g\rangle$를 **켓(ket)**이라고 부릅니다.

사실 브라와 켓은 따로 쓰면 각각 행벡터와 열벡터인데 양자역학에서는 켓은 계의 상태를 의미하며 브라는 (벡터를 스칼라로 보내는)선형 사상을 의미합니다. 위의 식에서 볼 수 있듯이 브라와 켓이 결합하면 내적이 됩니다. 더 깊이 들어가면 Hilbert 공간 같은 얘기를 해야 하는데 그냥 넘어가도록 하죠... 아무튼 지금은 함수의 내적을 저렇게 표기하면 된다 정도로 아시면 됩니다.

만약 중간에 연산자가 끼어 있는 경우에는 이렇게 씁니다.

$$
\langle f | \hat{A} | g \rangle \equiv \int f^* \hat{A} g \, dx
$$

그러면 에르미트 조건을 훨씬 짧고 깔끔하게 쓸 수 있습니다.

$$
\langle f | \hat{A} | g \rangle = \langle \hat{A} f | g \rangle
$$

연산자가 벽을 넘어 왼쪽으로 옮겨갈 수 있다는 뜻이죠. 앞으로는 이 브라켓 표기를 계속 쓰겠습니다.

0장에서 봤던 대응표 일부를 이 표기로 다시 써봅시다.

| 유한차원 | 무한차원 |
|---|---|
| 벡터 $\mathbf{v}$ | 켓 $| f \rangle$ |
| 내적 $\langle \mathbf{u}, \mathbf{v} \rangle$ | $\langle f | g \rangle$ |
| 행렬 $A$ | 연산자 $\hat{A}$ |
| 에르미트 행렬 | 에르미트 연산자 |
| $A\mathbf{v} = \lambda \mathbf{v}$ | $\hat{A} | f \rangle = a | f \rangle$ |


## 왜 운동량 연산자에는 허수가 붙나

이제 앞 절의 숙제를 풀어봅시다. 운동량 연산자에 $-i\hbar$가 붙은 진짜 이유입니다.

먼저 $i$ 없이 미분만 있다고 해봅시다. 즉 $\hat{A} = d/dx$인 경우죠. 에르미트인지 확인하려면 정의에 넣어보면 됩니다.

$$
\int f^* \left( \frac{d}{dx}g \right) dx
$$

부분적분 공식을 쓰면 다음과 같습니다.

$$
= \Big[f^*g\Big]_{-\infty}^{\infty} - \int \left( \frac{d}{dx} f^* \right) g \, dx
$$

파동함수는 무한대에서 0이 되어야 하니까 첫 항은 0이 되겠네요. (파동함수가 무한대에서 0이 되지 않으면 적분값이 발산하게 됩니다. 3.4절의 조건을 떠올려보세요.) 그러면 이렇게 됩니다.

$$
\int f^* \left( \frac{d}{dx}g \right) dx = -\int \left(\frac{d}{dx} f \right)^* g \, dx
$$

앗, 에르미트가 되려면 부호가 같아야 하는데 반대가 되어버렸습니다! 이런 경우를 반에르미트(anti-Hermitian)라고 부릅니다.

만약 처음에 $-i\hbar$가 곱해져 있었다면 어떻게 될까요? 부분적분을 하면서 부호가 반대가 되지만 켤레를 씌울 때 다시 부호가 바뀌니 제자리로 돌아옵니다.

$$
\int f^* \left(-i\hbar \frac{d}{dx} g \right) dx = \int \left(i\hbar \frac{d}{dx} f^* \right) g \, dx = \int \left(-i\hbar \frac{d}{dx} f \right)^* g \, dx
$$

$\hat{p} = -i\hbar \frac{d}{dx}$는 에르미트 연산자라는 것을 확인했습니다!

0.2.3절의 팁이 떠오르시나요? 그때 이렇게 썼었죠.

> $D_1$은 위쪽 대각선이 $+$, 아래쪽이 $-$인 반대칭 행렬입니다. 그래서 `eigh`를 쓸 수 없죠. 반대칭 행렬에 $i$를 곱하면 에르미트가 됩니다.

행렬에서 본 것과 정확히 같은 이야기입니다. 그때는 수치계산 요령처럼 보였는데, 이제 그것이 "운동량이 실수 값으로 측정되기 위한 조건"이었다는 게 드러났습니다.


## 코드로 확인해보기

물론 실제로 컴퓨터로 계산을 하려면 연산자를 행렬로 바꿔줘야 합니다. 0.2.3의 미분 행렬로 세 연산자를 만들고 에르미트 행렬인지 검사하겠습니다.

```python
import numpy as np
from scipy.constants import hbar, m_e, e

L, N = 1e-9, 1000
dx = L / (N + 1)
x = np.linspace(0, L, N + 2)[1:-1]

D1 = (np.eye(N, k=1) - np.eye(N, k=-1)) / (2 * dx)
D2 = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2

X = np.diag(x)  # 위치
P = -1j * hbar * D1  # 운동량
T = -(hbar**2) / (2 * m_e) * D2  # 운동에너지


def check(M, name):
    diff = np.max(np.abs(M - M.conj().T))
    print(
        f"{name:14s} |M - M^dagger| = {diff:.3e}   에르미트? {np.allclose(M, M.conj().T)}"
    )


check(X, "x")
check(P, "p")
check(T, "T")
check(D1, "D1 (without i)")
```
```
x              |M - M^dagger| = 0.000e+00   에르미트? True
p              |M - M^dagger| = 0.000e+00   에르미트? True
T              |M - M^dagger| = 0.000e+00   에르미트? True
D1 (without i) |M - M^dagger| = 1.001e+12   에르미트? False
```

위에서 본 것처럼 $D_1$은 에르미트가 아닙니다. 그런데 여기에 $-i\hbar$를 곱한 $\hat{p}$는 에르미트죠.

고유값도 봅시다.

```python
for M, name in [(X, "x"), (P, "p"), (T, "T"), (D1, "D1")]:
    ev = np.linalg.eigvals(M)
    print(
        f"{name:3s}: 실수부 최대 {np.max(np.abs(ev.real)):.3e}  허수부 최대 {np.max(np.abs(ev.imag)):.3e}"
    )
```
```
x  : 실수부 최대 9.990e-10  허수부 최대 0.000e+00
p  : 실수부 최대 1.056e-22  허수부 최대 0.000e+00
T  : 실수부 최대 2.447e-14  허수부 최대 0.000e+00
D1 : 실수부 최대 2.747e-04  허수부 최대 1.001e+12
```

$D_1$의 고유값은 순허수입니다. (실수부의 최대값이 허수부에 비해서 무시할 만한 정도입니다.) 그림으로 보면 더욱 분명해집니다.

```python
import matplotlib.pyplot as plt

fig, ax = plt.subplots(figsize=(6, 5))
evD1 = np.linalg.eigvals(D1)
evP = np.linalg.eigvals(P)

ax.plot(evD1.real / 1e10, evD1.imag / 1e10, "o", ms=4, color="crimson", label="D_1")
ax.plot(evP.real / 1e-25, evP.imag / 1e-25, "s", ms=4, color="steelblue", label="p")
ax.axhline(0, color="black", lw=0.8)
ax.axvline(0, color="black", lw=0.8)
ax.set_xlabel("Re (scaled)")
ax.set_ylabel("Im (scaled)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![연산자와 고유값](/assets/image-53.png)

$D_1$의 고유값은 세로축(허수축)에, $\hat{p}$의 고유값은 가로축(실수축)에 놓입니다. 복소평면에서 $-i\hbar$를 곱한다는 것은 90도 회전시키는 것이고, 그래서 허수축의 값들이 실수축으로 옮겨간 셈이죠.


## 고유함수는 직교한다

에르미트의 두 번째 선물입니다. 서로 다른 고유값에 속한 고유함수는 직교합니다. 상자 속 입자로 확인해봅시다.

```python
def psi(n):
    return np.sqrt(2 / L) * np.sin(n * np.pi * x / L)


print("      " + "".join(f"{'n='+str(n):>10}" for n in range(1, 4)))
for m in range(1, 4):
    row = "".join(f"{np.trapezoid(psi(m)*psi(n), x):10.5f}" for n in range(1, 4))
    print(f"m = {m} {row}")
```
```
             n=1       n=2       n=3
m = 1    1.00000   0.00000  -0.00000
m = 2    0.00000   1.00000   0.00000
m = 3   -0.00000   0.00000   1.00000
```

$$
\langle\psi_m|\psi_n\rangle = \delta_{mn}
$$

0.1.2에서 정규직교기저를 정의할 때 나왔던 그 식입니다. 그리고 0.1.5에서 사인 함수들이 직교한다는 것도 이미 확인했었죠. 그때는 순수하게 수학적인 성질이었는데, 이제 물리적인 의미가 생겼습니다. Hamiltonian이 에르미트이기 때문입니다.


## 다음 이야기

그런데 이게 왜 중요할까요? 0.1.5에서 직교하는 완전기저가 있으면 임의의 함수를 이 기저를 이용해서 펼칠 수 있었다고 했던 것 기억나시나요?

$$
\psi = \sum_n c_n\psi_n \qquad c_n = \langle \psi_n | \psi \rangle
$$

그런데 여기에서 나오는 계수 $c_n$들은 무엇을 의미하는 걸까요? 그리고 $\langle f | \hat{A} | g \rangle$ 같은 표기는 대체 무슨 의미를 갖는 걸까요?


## 확인 문제

1. 위치 연산자 $\hat{x}$가 에르미트임을 정의에 넣어 확인해보세요. (힌트: $x$는 실수입니다)
2. $\hat{A} = \frac{d^2}{dx^2}$은 에르미트일까요?
3. $\hat{p}$의 정의에서 $-i\hbar$ 대신 $+i\hbar$를 써도 에르미트일까요? 그렇다면 왜 부호를 $-$로 정했을까요?
4. 두 에르미트 연산자의 곱 $\hat{A}\hat{B}$는 에르미트일까요? $(\hat{A}\hat{B})^\dagger = \hat{B}^\dagger\hat{A}^\dagger$임을 이용해보세요.