# 5.4. 교환자와 불확정성

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/05-04.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞 절에서 운동량 연산자의 분산은 0이 아니라고 했습니다. 평균은 0이었지만 제곱한 값의 평균이 0이 아니었거든요. 그런데 코드로 계산해보기는 했지만 구체적으로 이게 무슨 의미인지는 잘 와닿지가 않습니다... 연산자의 분산의 크기가 대체 왜 중요할까요? 사실 여기에는 더 심오한 것이 숨어 있습니다.


## 연산자와 표준편차

연산자의 분산을 계산한 것과 마찬가지로 연산자의 **표준편차(standard deviation)**도 계산할 수 있습니다!

$$
\sigma_A = \sqrt{\langle A^2\rangle - \langle A\rangle^2}
$$

표준편차는 통계학에서 데이터가 평균에서 얼마나 벗어나 있는지를 재는 양입니다. 2.4절에서 파동 묶음의 폭을 잴 때 이미 쓴 적이 있죠. $\sigma_A = 0$이면 모든 데이터가 평균에 있다는 뜻입니다. 즉, 측정할 때마다 매번 같은 값이 나온다는 것이고, 앞 절에서 봤듯이 그건 계가 $\hat{A}$의 고유상태에 있다는 뜻입니다.

상자 속 입자에서 위치 연산자와 운동량 연산자의 표준편차를 계산해봅시다.

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import hbar, m_e, e

L, N = 1e-9, 3000
dx = L / (N + 1)
xg = np.linspace(0, L, N + 2)[1:-1]
D1 = (np.eye(N, k=1) - np.eye(N, k=-1)) / (2 * dx)
D2 = (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2


def psi(n):
    return np.sqrt(2 / L) * np.sin(n * np.pi * xg / L)


print(f"{'n':>2} {'σx/L':>9} {'σp (hbar/L)':>13} {'σx * σp (hbar)':>15}")
for n in [1, 2, 3, 4, 5]:
    p = psi(n)
    x_mean = np.trapezoid(p * xg * p, xg)
    x2_mean = np.trapezoid(p * xg**2 * p, xg)
    p2_mean = np.trapezoid(p * (-(hbar**2) * (D2 @ p)), xg)

    sigma_x = np.sqrt(x2_mean - x_mean**2)
    sigma_p = np.sqrt(p2_mean)  # <p> = 0 이므로
    print(f"{n:2d} {sigma_x/L:9.4f} {sigma_p/(hbar/L):13.4f} {sigma_x * sigma_p/hbar:15.4f}")
```
```
 n      σx/L   σp (hbar/L)  σx * σp (hbar)
 1    0.1808        3.1416          0.5679
 2    0.2658        6.2832          1.6703
 3    0.2788        9.4248          2.6272
 4    0.2831       12.5664          3.5580
 5    0.2851       15.7079          4.4790
```

$\sigma_x \sigma_p$가 항상 $0.5\hbar$보다 큽니다. 그리고 바닥상태가 $0.568\hbar$로 가장 가깝죠. 2.4절에서 봤던 불확정성 원리와 같은 결론입니다.

그림으로 보면 더 분명합니다.

```python
ns = np.arange(1, 9)
Sigma_xs, Sigma_ps, prods = [], [], []

for n in ns:
    p = psi(n)
    x_mean = np.trapezoid(p * xg * p, xg)
    x2_mean = np.trapezoid(p * xg**2 * p, xg)
    p2_mean = np.trapezoid(p * (-(hbar**2) * (D2 @ p)), xg)
    sigma_x, sigma_p = np.sqrt(x2_mean - x_mean**2), np.sqrt(p2_mean)
    Sigma_xs.append(sigma_x / L)
    Sigma_ps.append(sigma_p / (hbar / L))
    prods.append(sigma_x * sigma_p / hbar)

fig, (a1, a2) = plt.subplots(1, 2, figsize=(11, 5))

a1.plot(ns, Sigma_xs, "o-", color="crimson", lw=2)
a1.set_xlabel("n")
a1.set_ylabel("σx / L", color="crimson")
a1.tick_params(axis="y", labelcolor="crimson")
a1.set_ylim(0, 0.32)
a1.grid(alpha=0.3)
a1b = a1.twinx()
a1b.plot(ns, Sigma_ps, "s--", color="steelblue", lw=2)
a1b.set_ylabel("σp / (hbar/L)", color="steelblue")
a1b.tick_params(axis="y", labelcolor="steelblue")

a2.plot(ns, prods, "o-", color="darkorange", lw=2, label="σx σp")
a2.axhline(0.5, color="black", ls="--", lw=1.6, label="lower bound hbar/2")
a2.set_xlabel("n")
a2.set_ylabel("σx σp / hbar")
a2.legend(fontsize=9)
a2.grid(alpha=0.3)

plt.tight_layout()
plt.show()
```
![연산자의 표준편차](/assets/image-56.png)

$\sigma_x$는 금방 포화됩니다. $n$이 커져도 0.29$L$ 근처에서 멈추죠. 당연합니다. 전자가 상자 안에 갇혀 있으니 위치의 흔들림이 상자 크기를 넘을 수 없거든요. 반면에 $\sigma_p$는 $n$에 비례해 계속 커집니다. $\langle p^2\rangle = (\frac{n\pi\hbar}{L})^2$이었으니까요. 그러니 둘의 곱은 커질 수밖에 없습니다.

*아니, 그건 당연해 보이는데 그래서 뭐가 중요하다는 거예요?*


## 위치와 운동량의 교환자

이제 위치 연산자와 운동량 연산자의 교환자를 계산해봅시다. 수치적으로 말고 손으로 직접이요! 임의의 함수 $\psi$에 작용시켜보면 되죠.

$$
\hat{x}\hat{p}\,\psi = x\left(-i\hbar\frac{d\psi}{dx}\right) = -i\hbar\,x\frac{d\psi}{dx}
$$

$$
\hat{p}\hat{x}\,\psi = -i\hbar\frac{d}{dx}\left(x\psi\right) = -i\hbar\left(\psi + x\frac{d\psi}{dx}\right)
$$

곱의 미분법 때문에 아래쪽에 항이 하나 더 생겼습니다. 따라서 둘의 차이는 이렇게 됩니다.

$$
[\hat{x},\hat{p}]\,\psi = \hat{x}\hat{p}\psi - \hat{p}\hat{x}\psi = i\hbar\,\psi
$$

이 식은 $\psi$가 무엇이든 상관없이 성립합니다. 그러니 교환자 형식으로 쓰면 이렇게 되네요.

$$
[\hat{x},\hat{p}] = i\hbar
$$

코드로도 확인해봅시다.

```python
X = np.diag(xg)  # 위치 연산자 (대각행렬)
P = -1j * hbar * D1  # 운동량 연산자

C = X @ P - P @ X  # 교환자

for n in [1, 2, 3, 4]:
    p = psi(n)
    value = np.trapezoid(p * (C @ p), xg)
    print(f"n = {n}: <[x, p]> = {value / (1j * hbar):.4f} x (i hbar)")
```
```
n = 1: <[x, p]> = 1.0000+0.0000j x (i hbar)
n = 2: <[x, p]> = 1.0000+0.0000j x (i hbar)
n = 3: <[x, p]> = 1.0000+0.0000j x (i hbar)
n = 4: <[x, p]> = 1.0000+0.0000j x (i hbar)
```

모든 상태에서 위치 연산자와 운동량 연산자의 교환자는 항상 $i\hbar$입니다.


## 교환하는 연산자, 교환하지 않는 연산자

모든 연산자가 교환하지 않는 것은 아닙니다. 몇 가지를 확인해봅시다.

```python
H = -(hbar**2) / (2 * m_e) * D2
V = np.diag(xg**2)


def commutator(A, B):
    return A @ B - B @ A


for A, B, name in [(X, V, "[x, V(x)]"), (H, H, "[H, H]  "), (X, H, "[x, H]  ")]:
    print(f"{name}: 최대 크기 = {np.max(np.abs(commutator(A, B))):.3e}")
```
```
[x, V(x)]: 최대 크기 = 0.000e+00
[H, H]  : 최대 크기 = 0.000e+00
[x, H]  : 최대 크기 = 1.832e-26
```

$\hat{x}$와 $\hat{V}(x)$는 교환합니다. 둘 다 곱셈 연산자니 순서를 바꿔도 같은 수를 곱하는 것뿐이죠. 자기 자신과도 당연히 교환하고요. 반면 $\hat{x}$와 $\hat{H}$는 교환하지 않습니다. $\hat{H}$ 안에 미분이 들어 있으니까요.

여기서 중요한 사실이 나옵니다.

> 두 연산자가 교환하면 공통의 고유함수를 가질 수 있고, 두 물리량을 동시에 확정된 값으로 알 수 있다.

왜 그럴까요? $\hat{A}\psi = a\psi$이고 $\hat{B}\psi = b\psi$인 상태가 있다면

$$
\hat{A}\hat{B}\psi = \hat{A}(b\psi) = ab\psi, \qquad \hat{B}\hat{A}\psi = \hat{B}(a\psi) = ba\psi
$$

$a$와 $b$는 그냥 숫자니 $ab = ba$이고, 두 결과가 같습니다. 즉 공통의 고유함수가 존재하면 교환자는 그 상태에서 0이죠. 이 성질은 앞으로도 계속 쓰입니다.

반대로 교환자가 0이 아니면 공통 고유함수가 있을 수 없습니다. $[\hat{x},\hat{p}] = i\hbar \neq 0$이니 위치와 운동량을 동시에 확정할 수 있는 상태는 존재하지 않습니다. 불확정성 원리가 여기서 나오는 것입니다!


## 불확정성 원리와 교환자

일반적으로 연산자의 표준편차와 교환자 사이에는 다음 부등식이 성립합니다. 자세한 유도 과정은 생략하도록 하겠습니다.

$$
\sigma_A^2 \sigma_B^2 \ge \left( \frac{1}{2i} \left\langle[\hat{A},\hat{B}]\right\rangle \right)^2
$$

교환자가 두 물리량의 표준편차의 곱에 하한을 준다는 것이 핵심입니다. 교환자가 크면 하한도 높고, 0이면 하한이 없죠.

위치와 운동량을 넣어봅시다. $[\hat{x},\hat{p}] = i\hbar$으로 일정하니까 기댓값도 상수고요.

$$
\sigma_x^2 \sigma_p^2 \ge \left( \frac{1}{2i} (i\hbar) \right)^2 = \frac{\hbar^2}{4} \qquad \Longrightarrow \qquad \sigma_x \sigma_p \ge \frac{\hbar}{2}
$$

Heisenberg의 불확정성 원리가 나왔습니다! 2.4절에서는 여러 파장을 섞어 묶음을 만들고 그 폭을 재서 얻었죠. 지금은 연산자의 순서를 바꿔본 것뿐인데 같은 곳에 도달했습니다. 완전히 다른 두 길이 한 지점에서 만난 겁니다.

여기서 맨 위 코드의 계산 결과를 잠깐 다시 보시죠.

```
n = 1: σx * σp = 0.5679 hbar
```

확실히 $0.5\hbar$보다 크긴 한데, 왜 정확히 $0.5\hbar$가 아닐까요?

2.4절에서 가우스 함수 형태의 파동을 섞었을 때는 정확히 $\hbar/2$가 나왔습니다. 폭이 얼마든 상관없이요. 즉, $\hbar/2$는 가우스 함수 모양일 때만 달성되는 값입니다. 사인 함수는 상자 벽에서 뚝 끊기는 모양이라 조금 손해를 보죠. 그래서 조금 더 큰 값이 나온 겁니다.

이 성질 때문에 가우스 함수를 **최소 불확정성 상태(minimum uncertainty state)**라고 부릅니다. 그리고 이것이 나중에 중요해집니다. 6장의 조화진동자 바닥상태가 정확히 가우스 모양이거든요. 16장에서 계산화학이 기저함수로 가우스 함수를 쓰는 것도 무관하지 않고요.


## 정리하며

조금 추상적으로 보일 수도 있었던 5장이 끝났습니다. 다음 장으로 넘어가기 전에 양자역학의 중요한 기본 가정들을 정리하고 가겠습니다.  

**가정 1**: 양자역학적 계의 상태는 파동함수에 의해 완전히 결정되며 파동함수의 제곱은 입자가 발견될 확률밀도를 나타낸다.  
**가정 2**: 측정 가능한 모든 물리량은 어떤 에르미트 연산자와 대응되고, 측정값은 항상 고유값 중 하나로 나타난다.  
**가정 3**: 연산자 $A$에 대응하는 물리량을 반복해서 측정했을 때의 평균값은 $\langle A \rangle$로 주어진다.  
**가정 4**: 계의 파동함수의 시간에 따른 변화는 시간 의존 Schrödinger 방정식을 따른다.  

마지막 가정은 지금 당장 중요하지는 않지만 나중에 분자 동역학이나 분자가 빛을 방출하거나 흡수하는 현상을 설명할 때 시간에 의존하는 Schrödinger 방정식을 써야 할 때가 있다는 정도만 알아두시면 됩니다.  

다음 장부터는 다시 실제 문제로 돌아가서 계산을 해보도록 하죠.


## 확인 문제

1. $[\hat{p}, \hat{x}]$는 얼마일까요? $[\hat{x},\hat{p}]$와 부호가 어떻게 다른가요?
2. $[\hat{x}, \hat{x}^2]$을 계산해보세요. 곱셈 연산자끼리는 왜 항상 교환할까요?
3. $\hat{T}$와 $\hat{V}$는 교환할까요? 코드로 확인해보고, 그 결과가 $[\hat{x},\hat{H}] \neq 0$과 어떻게 연결되는지 생각해보세요.
4. 상자 속 입자의 $n$을 100까지 늘려보세요. $\sigma_x \sigma_p$가 어떻게 변하나요? 1.4절의 대응 원리와 관련지어 설명해보세요.