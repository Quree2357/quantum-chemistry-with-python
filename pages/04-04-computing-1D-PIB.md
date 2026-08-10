# 4.4. 컴퓨터로 풀어보기

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/04-04.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

지금까지 상자 속 입자 모델을 손으로 풀어봤습니다. 그런데 사실 모두 상자 안에서 퍼텐셜이 0으로 일정해서 손으로 풀어도 그렇게 어렵지 않았습니다. 퍼텐셜이 조금만 복잡해져도 손으로는 답이 안 나옵니다. 이제 컴퓨터에게 넘길 차례입니다.

그리고 좋은 소식이 있습니다. 여러분이 0장에서부터 가지고 다녔던 무기들이 다시 활약합니다.


## 이도류 스킬 하나면 된다

3.4절의 마지막 부분에서 우리는 Schrödinger 방정식이 고유값 문제인 것을 확인했고 어떻게 푸는지도 알고 있다고 했습니다.

$$
\hat{H} = -\frac{\hbar^2}{2m}\frac{d^2}{dx^2} + V(x)
\qquad\Longrightarrow\qquad
\mathbf{H} = -\frac{\hbar^2}{2m}D_2 + \mathrm{diag}(V)
$$

0.2.3에서 만든 2차 미분 행렬 $D_2$의 대각선에 퍼텐셜을 얹으면 Hamiltonian 행렬이 됩니다. 퍼텐셜이 대각행렬인 이유는 곱셈 연산자이기 때문입니다. $V(x)\psi(x)$는 각 격자점에서 그 자리의 값을 곱하는 것뿐이니, 다른 격자점과 섞이지 않죠.

코드로는 이렇게 계산할 수 있습니다.

```python
import numpy as np
from scipy.constants import hbar, m_e, e, h


def second_derivative_matrix(N, dx):
    """0.2.3절에서 만든 그 행렬"""
    return (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2


def solve(N, L, potential=None):
    """길이 L인 구간을 N개 격자로 나눠 Schrödinger 방정식을 푼다"""
    dx = L / (N + 1)
    x = np.linspace(0, L, N + 2)[1:-1]  # 양 끝에서는 psi=0 이므로 제외
    V = np.zeros(N) if potential is None else potential(x)

    H = -(hbar**2) / (2 * m_e) * second_derivative_matrix(N, dx) + np.diag(V)
    E, psi = np.linalg.eigh(H)

    return x, E / e, psi
```

핵심은 함수 안의 두 줄입니다. 행렬을 조립하고 `eigh`를 부르는 것, 그게 전부죠. 검으로 두 번 베는 이도류 스킬을 썼더니 그냥 끝나버렸습니다.

$$
\underbrace{-\frac{\hbar^2}{2m}D_2 + \mathrm{diag}(V)}_{\texttt{H = ...}}
\qquad\longrightarrow\qquad
\underbrace{\hat{H}\psi = E\psi}_{\texttt{E, psi = eigh(H)}}
$$

0.1.3에서 `eigh`를 배울 때 "이 책의 거의 모든 문제는 결국 `eigh` 한 줄로 끝난다"고 했었죠. 그 순간이 왔습니다.


## 상자 속 입자 다시 풀어보기

먼저 답을 아는 문제부터 풀어봐야 합니다. 코드를 믿기 전에 검증해야 하니까요. 4.1절에서 상자 속 입자에 대해 계산했던 결과는 이랬습니다.

$$
E_n = \frac{n^2h^2}{8mL^2}
$$

$L=1$ nm 상자에 이 코드를 돌려봅시다.

```python
L = 1e-9


def exact(n, L=L):
    return n**2 * h**2 / (8 * m_e * L**2) / e


print(f"{'N':>7} " + "".join(f"{'n='+str(n):>11}" for n in range(1, 5)))
for N in [10, 50, 200, 1000]:
    x, E, psi = solve(N, L)
    print(f"{N:7d} " + "".join(f"{E[n-1]:11.5f}" for n in range(1, 5)))
print(f"{'해석해':>5} " + "".join(f"{exact(n):11.5f}" for n in range(1, 5)))
```
```
      N        n=1        n=2        n=3        n=4
     10    0.37348    1.46367    3.18224    5.38997
     50    0.37591    1.50222    3.37465    5.98610
    200    0.37602    1.50400    3.38365    6.01452
   1000    0.37603    1.50412    3.38425    6.01640
  해석해    0.37603    1.50412    3.38427    6.01648
```

답이 맞습니다! 격자를 촘촘히 할수록 해석해에 가까워집니다. $N=1000$이면 소수점 넷째 자리까지 일치하네요.

여기서 잠깐 멈춰서 무슨 일이 일어났는지 봅시다. 우리는 컴퓨터에게 사인 함수도, 양자수도, $n^2$ 규칙도 알려주지 않았습니다. 격자와 퍼텐셜만 주고 행렬의 고유값을 구하라고 했을 뿐이죠. 그런데 4.1절에서 손으로 유도했던 결과가 똑같이 나왔습니다.

파동함수도 확인해봅시다.

```python
import matplotlib.pyplot as plt

x, E, psi = solve(200, L)

fig, ax = plt.subplots(figsize=(6, 4))
for n, col in zip([1, 2, 3], ["red", "orange", "green"]):
    v = psi[:, n - 1]
    if v[len(v) // 4] < 0:  # 부호 관례를 맞춰준다
        v = -v
    ax.plot(x * 1e9, v / np.abs(v).max(), lw=3, color=col, label=f"n={n}")
    ax.plot(
        x * 1e9, np.sin(n * np.pi * x / L), ls="--", lw=2, color="black", alpha=1
    )

ax.set_xlabel("x (nm)")
ax.hlines(0, 0, L * 1e9, color="black", lw=1)
ax.vlines([0, L * 1e9], -1.2, 1.2, color="black", lw=1)
ax.set_ylim(-1.2, 1.2)
ax.set_yticks([])
ax.legend(fontsize=9)
plt.show()
```
![컴퓨터로 1-D PIB 풀기](/assets/image-47.png)

검은 점선이 해석해 $\sin(n\pi x/L)$입니다. 색깔 곡선과 완전히 겹쳐 있죠.

[[TIP]]
`eigh`가 돌려주는 고유벡터의 부호는 정해져 있지 않습니다. $\psi$가 해라면 $-\psi$도 해거든요. 그래서 코드에서 부호를 뒤집는 줄을 넣었습니다.
그리고 고유벡터는 벡터로서 크기가 1로 규격화되어 있는데, 이건 $\int|\psi|^2dx=1$과는 다릅니다. 격자 간격 $dx$를 곱해야 하죠. 정확히는 $\psi_{\text{물리}} = \psi_{\text{eigh}}/\sqrt{dx}$입니다. 확률을 계산할 때는 이 점을 챙겨야 합니다.
[[/TIP]]


## 얼마나 촘촘해야 하는가

격자를 늘리면 좋아진다는 건 알겠는데, 얼마나 늘려야 할까요? 그리고 준위마다 다를까요?

```python
Ns = [5, 10, 20, 50, 100, 200, 500, 1000]

fig, ax = plt.subplots(figsize=(6.5, 4.5))
for n, col in zip([1, 2, 3, 4], ["red", "orange", "green", "blue"]):
    errors = [abs(solve(N, L)[1][n - 1] - exact(n)) / exact(n) * 100 for N in Ns]
    ax.loglog(Ns, errors, "o-", color=col, label=f"n={n}")

ax.set_xlabel("N (grid points)")
ax.set_ylabel("relative error (%)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3, which="both")
plt.show()
```
![격자 크기에 따른 오차](/assets/image-48.png)

네 직선이 나란히 아래로 내려갑니다. 기울기를 보면 $-2$인데, 오차가 $N^{-2}$에 비례한다는 뜻이죠. 0.2.1에서 중심차분이 $O(h^2)$였으니 당연한 결과입니다.

그런데 네 직선의 높이가 다릅니다. 숫자로 보면 이렇습니다.

```python
for N in [10, 50, 200]:
    x, E, psi = solve(N, L)
    line = "  ".join(
        f"n={n}: {abs(E[n-1]-exact(n))/exact(n)*100:7.4f}%" for n in range(1, 5)
    )
    print(f"N={N:5d} : {line}")
```
```
N=   10 : n=1:  0.6779%  n=2:  2.6895%  n=3:  5.9698%  n=4: 10.4134%
N=   50 : n=1:  0.0316%  n=2:  0.1264%  n=3:  0.2843%  n=4:  0.5049%
N=  200 : n=1:  0.0020%  n=2:  0.0081%  n=3:  0.0183%  n=4:  0.0326%
```

위쪽 준위일수록 오차가 큽니다. $N=10$일 때 바닥상태는 0.7%인데 $n=4$는 10%가 넘죠.  
이유는 간단합니다. $n$이 클수록 파동함수의 출렁임이 촘촘해지는데, 격자가 그걸 따라가지 못하는 겁니다. 마디가 $n-1$개인 함수를 표현하려면 최소한 그보다 훨씬 촘촘한 격자가 필요하죠.

> 원칙: 격자 계산에서는 아래쪽 준위 몇 개만 믿는다.

바닥상태와 그 위 몇 개는 잘 맞지만, 위로 갈수록 의심해야 합니다. 이건 격자 방법만의 문제가 아닙니다. 유한한 자원으로 무한차원 문제를 푸는 모든 방법이 같은 성질을 가집니다. 나중에도 계속 같은 이야기가 반복될 겁니다. 격자가 촘촘해지면 정확도는 높아지지만 계산량도 같이 늘어나거든요.


## 어떤 퍼텐셜 공격이 들어와도 막을 수 있다

검증이 끝났으니 본론입니다. 이 코드에서 $V$만 갈아끼우면 다른 문제가 됩니다.  
몇 가지 퍼텐셜을 넣어봅시다. 계단, 기울어진 바닥, 그리고 살짝 튀어나온 턱입니다.

```python
L, N = 2e-9, 500

potentials = [
    ("flat box", lambda x: np.zeros_like(x)),
    ("step", lambda x: np.where(x > L / 2, 1.0 * e, 0.0)),
    ("triangular", lambda x: 2.0 * e * x / L),
    ("barrier", lambda x: np.where(np.abs(x - L / 2) < L / 10, 1.5 * e, 0.0)),
]

fig, axes = plt.subplots(1, 4, figsize=(14, 4))

for ax, (title, Vf) in zip(axes, potentials):
    x, E, psi = solve(N, L, Vf)
    ax.plot(x * 1e9, Vf(x) / e, color="black", lw=1.6)  # 퍼텐셜
    for k in range(4):
        v = psi[:, k]
        if v[np.argmax(np.abs(v))] < 0:
            v = -v
        ax.hlines(E[k], 0, L * 1e9, color="gray", lw=0.7, ls=":")
        ax.plot(x * 1e9, E[k] + 0.9 * v / np.abs(v).max(), lw=1.6)
    ax.set_xlim(0, L * 1e9)
    ax.set_ylim(-1, E[4] + 1.5)
    ax.set_xlabel("x (nm)")
    ax.set_title(title, fontsize=11)

axes[0].set_ylabel("E (eV)")
plt.show()
```
![퍼텐셜에 따른 모드](/assets/image-49.png)

각 준위의 에너지 높이에 파동함수를 얹어 그렸습니다. 교재에서 흔히 보는 그림이죠. 바뀐 것은 퍼텐셜 하나뿐입니다. 그런데 결과가 전혀 다릅니다.  
퍼텐셜이 높은 곳에서는 확실히 전자가 발견될 확률이 작아 보입니다. 이것은 모든 그림에서 똑같이 나타나죠. 그런데 마지막 그림이 조금 이상합니다.  
전자의 에너지가 언덕의 높이보다 낮은 상태에 있는데도 파동함수가 언덕 안에서 0이 아닙니다! 고전적인 접근에서는 전자는 언덕을 넘어갈 수 없습니다. 무언가 알 수 없는 일이 벌어진 것 같네요. 다음 절에서 이것에 대해 이야기해 보죠.

[[TIP]]
왜 3차원 문제는 코드로 풀어보지 않았을까요? 바로 격자의 크기 때문입니다. 3차원에서 각 축을 $N$개로 나누면 격자점이 $N^3$개가 되고, Hamiltonian은 $N^3 \times N^3$ 행렬이 됩니다. 한 축당 50개씩만 잡아도 $125{,}000 \times 125{,}000$이라 `eigh`로는 감당이 안 되죠. 메모리만 100 GB가 넘습니다.
다행히 4.3절의 3차원 상자의 경우에는 변수분리로 깔끔하게 풀렸습니다. 하지만 변수분리가 안 되는 3차원 문제라면 다른 방법이 필요하죠. 실제 분자를 다룰 때는 그런 방법들을 쓰는데, 16장에서 그 이야기를 하겠습니다.
[[/TIP]]


## 확인 문제

1. `solve` 함수로 $L=2$ nm 상자를 풀어보세요. 바닥상태 에너지가 $L=1$ nm일 때의 몇 분의 1인가요? 4.1절에서 배운 관계와 맞나요?
2. 계단 퍼텐셜에서 계단의 높이를 5 eV로 올려보세요. 아래쪽 준위의 파동함수가 어떻게 달라지나요?
3. 가운데 언덕의 높이를 1.5 eV에서 5 eV로 올려보세요. 언덕 안으로 스며드는 정도가 어떻게 달라지나요?
4. 위쪽 준위의 오차가 큰 이유를 설명해보세요. $n=10$인 상태를 제대로 계산하려면 격자를 몇 개쯤 잡아야 할까요?