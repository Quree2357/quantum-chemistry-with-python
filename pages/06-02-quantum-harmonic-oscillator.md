# 6.2. 양자 조화 진동자

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/06-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

양자역학적 조화 진동자는 앞에서 봤던 고전 조화 진동자 모델을 Schrödinger 방정식에 적용한 것입니다. 당연히 고전 모델보다는 풀기가 조금 더 까다롭지만 우리는 직접 풀어보지는 않고 결과에 조금 더 집중하기 위해 '날먹'을 해보도록 하죠.


## 용수철에 양자를 한 스푼

모든 조건을 고전 조화 진동자에서 그대로 가져옵니다. 퍼텐셜 에너지 $V(x)$는 $\frac{1}{2}kx^2$으로 표현되었었죠. Schrödinger 방정식에 대입해봅시다.

$$
\hat{H}\psi = - \frac{\hbar^2}{2m} \frac{d^2\psi}{dx^2} + \frac{1}{2}kx^2 \psi = E\psi
$$

4장에서 봤던 상자 속 입자 모델과 비교하면 두 가지 다른 점이 있습니다. 한 가지는 퍼텐셜 에너지 항이고, 다른 한 가지는 경계조건이 없다는 것입니다. 상자에는 가장자리에 벽이 있었지만 여기서는 아무것도 없습니다. 대신 3.4절의 조건이 그 역할을 해줍니다. 바로 파동함수의 규격화 조건입니다.  

식만 봐도 풀기 어려워 보입니다... 하지만 위에서 말했듯이 우리가 직접 풀지 않고 결과만 낼름 보겠습니다. 양자 조화 진동자는 Schrödinger 방정식을 해석적으로 풀 수 있는 몇 안 되는 경우여서 해가 정확히 알려져 있습니다.

$$
\psi_v(x) = A_v H_v \left( \sqrt{\alpha} x \right) e^{-\frac{\alpha x^2}{2}} \qquad \alpha = \frac{m\omega}{\hbar}
$$

$A_v$은 그냥 규격화 상수니까 별로 중요하지는 않습니다. 그 오른쪽에 보이는 $H_v$이라는 함수가 신경쓰이실 텐데요, 이것은 Hermite 다항식이라고 불리는 특수함수입니다. 처음 몇 개의 식만 볼까요?

$$
H_0(x) = 1, \quad H_1(x) = 2x, \quad H_2(x) = 4x^2 - 2, \quad H_3(x) = 8x^3 - 12x, \quad H_4(x) = 16x^4 - 48x^2 + 12, \cdots
$$

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-2.5, 2.5, 500)

H0 = 0 * x + 1
H1 = 2 * x
H2 = 4 * x**2 - 2
H3 = 8 * x**3 - 12 * x
H4 = 16 * x**4 - 48 * x**2 + 12
H5 = 32 * x**5 - 160 * x**3 + 120 * x

fig, axes = plt.subplots(figsize=(11, 7))

axes.plot(x, H0, lw=3, label='H_0')
axes.plot(x, H1, lw=3, label='H_1')
axes.plot(x, H2, lw=3, label='H_2')
axes.plot(x, H3, lw=3, label='H_3')
axes.plot(x, H4, lw=3, label='H_4')
axes.plot(x, H5, lw=3, label='H_5')
axes.axhline(0, lw=2, color='black')
axes.axvline(0, lw=2, ls='--', color='gray')
axes.set_xlim(-2.2, 2.2)
axes.set_ylim(-50, 50)
axes.legend(fontsize=9, loc='upper right')
axes.grid(alpha=0.5)

plt.show()
```
![Hermite 다항식](/assets/image-60.png)

n번째 Hermite 다항식은 n차 다항식이고 그림에서 볼 수 있듯이 n개의 서로 다른 근을 가집니다. 또, 홀수번째 Hermite 다항식은 기함수(odd function)이고 짝수번째 Hermite 다항식은 우함수(even function)네요.

[[TIP]]
대부분의 양자역학 교재는 양자 조화 진동자 문제를 **사다리 연산자(ladder operator)**를 사용해서 풉니다. 급수해 방법을 써서 직접 미분방정식을 풀 수도 있지만 특수한 연산자를 정의하고 그 관계들만으로 해를 얻어내는 방법이죠. Dirac이 발전시킨 것으로 알려져 있습니다. 이 방법은 미분방정식을 직접 풀지 않아도 될 뿐만 아니라 아름다운 구조도 있습니다. 정의한 연산자를 한 번 적용할 때마다 에너지 준위가 한 단계씩 오르거나 내려가서 모든 상태를 만들어낼 수 있거든요. 다만 이 책에서는 양자 조화 진동자의 해법에 대해서는 생략하겠습니다. 5.4절에서 배운 교환자가 실제로 위력을 발휘하는 자리이니 궁금하신 분들은 McQuarrie 5장이나 Griffiths 2장을 참고하세요.
[[/TIP]]


## 에너지는 어떻게 되나

위에서 파동함수가 어떻게 나오는지는 봤는데 에너지가 어떻게 나오는지는 안 봤군요! 양자 조화 진동자의 에너지는 다음과 같습니다.

$$
E_v = \hbar\omega \left( v + \frac{1}{2} \right)
$$

아니 파동함수는 그렇게 복잡했는데 에너지는 엄청 간단하게 나오네요? 상자 속 입자의 경우와 전혀 다른 모양입니다. 한번 비교해볼까요?
```python
fig, axes = plt.subplots(1, 2, figsize=(10, 5))

for n in range(1, 7):
    axes[0].hlines(n**2, 0, 1, lw=3, color="steelblue")
    axes[0].text(1.05, n**2, f"n={n}", va="center", fontsize=9)
for v in range(6):
    axes[1].hlines(v + 0.5, 0, 1, lw=3, color="crimson")
    axes[1].text(1.05, v + 0.5, f"v={v}", va="center", fontsize=9)

axes[0].set_xlim(0, 1.3)
axes[0].set_ylim(0, 40)
axes[0].set_xticks([])
axes[0].set_title("box: E ~ $n^2$", fontsize=11)
axes[0].set_ylabel("E (relative)")
axes[1].set_xlim(0, 1.3)
axes[1].set_ylim(0, 6)
axes[1].set_xticks([])
axes[1].set_title("oscillator: E ~ (v + 1/2)", fontsize=11)
axes[1].set_ylabel("E / ($\\hbar \\omega$) (relative)")
plt.show()
```
![상자 속 입자와 조화 진동자의 에너지 준위](/assets/image-61.png)

상자 속 입자의 경우에는 위로 갈수록 에너지 차이가 커지는데 조화 진동자의 경우에는 차이가 정확히 $\hbar\omega$로 일정합니다. 그리고 놀랍게도 조화 진동자에서는 축퇴가 없고, 모든 고유값은 각각 단 하나의 고유함수만 대응됩니다.


## 영점 에너지

상자 속 입자와 비교해서 더 독특한 차이점은 바로 **영점 에너지(zero-point energy)**의 존재입니다. 에너지 식에 $v=0$을 넣어도 에너지가 0이 되지 않죠.

$$
E_0 = \frac{1}{2}\hbar\omega
$$

분명히 퍼텐셜 에너지는 $x=0$일 때 0인데 왜 여기까지 내려가지 못하는 걸까요? 상자 속 입자에서도 영점 에너지가 존재하긴 했지만 이번에는 다른 이유입니다.  

상자 속 입자에서는 $n=0$이면 파동함수 자체가 0이 되어 전자가 사라져버리는 것이 문제였습니다. 하지만 이번에는 $v=0$이어도 파동함수가 0이 되지는 않죠. 이 경우에 문제가 되는 것은 불확정성 원리입니다. 전체 에너지가 0이 되려면 입자가 정확히 퍼텐셜 에너지 곡선의 최소점에 정지해 있어야 합니다. 하지만 불확정성 원리에 따르면 입자는 절대 정지해 있을 수 없습니다. 그렇지 않으면 $\Delta x$와 $\Delta p$가 모두 0이 되어버리니까요.  

이 영점 에너지는 실제로 어느 정도나 될까요? 수소 분자의 경우를 계산해봅시다. (여기서 질량 자리에는 원자 하나의 질량이 아니라 환산질량이 들어가게 됩니다! 왜 이렇게 되는지는 다음 절에서 살펴보겠습니다.)
```python
from scipy.constants import hbar, u, e

k = 572  # 수소 분자의 힘 상수 (N/m)
mu = 0.504 * u

omega = np.sqrt(k / mu)

print(f"omega         = {omega:.4e} rad/s")
print(f"hbar * omega  = {hbar * omega / e:.4f} eV")
print(f"영점 에너지    = {0.5 * hbar * omega / e:.4f} eV")
```
```
omega         = 8.2672e+14 rad/s
hbar * omega  = 0.5442 eV
영점 에너지    = 0.2721 eV
```
0.27 eV 정도의 영점 에너지를 갖네요. 온도를 절대영도까지 내려도 수소 분자는 이만큼의 진동 에너지를 갖게 됩니다. 작아 보이지만 화학에서는 절대 무시할 수 없는 크기입니다. 실제로 계산화학에서 반응 에너지를 구할 때 영점 에너지를 반드시 보정해줍니다.


## 이번에도 이도류가 간다

4.4절에서 만들었던 이도류 스킬(코드)을 또 써봅시다. 퍼텐셜 에너지만 바꾸면 됩니다.
```python
N, L = 2000, 5e-11
dx = 2 * L / (N + 1)
x = np.linspace(-L, L, N + 2)[1:-1]


def D2(N, dx):
    return (-2 * np.eye(N) + np.eye(N, k=1) + np.eye(N, k=-1)) / dx**2


def solve(N, L, potential=None):
    V = np.zeros(N) if potential is None else potential(x)
    H = -(hbar**2) / (2 * mu) * D2(N, dx) + np.diag(V)
    E, psi = np.linalg.eigh(H)

    return x, E / e, psi


V_harmonic = lambda x: 0.5 * k * x**2

x, E, psi = solve(N, L, V_harmonic)

print(f"{'v':>2} {'수치해 (eV)':>10} {'해석해 (eV)':>10} {'오차':>8}")
for v in range(6):
    exact = (v + 0.5) * hbar * omega / e
    print(f"{v:2d} {E[v]:12.5f} {exact:12.5f} {abs(E[v]-exact)/exact*100:9.4f}%")
```
```
 v   수치해 (eV)   해석해 (eV)       오차
 0      0.27208      0.27208    0.0000%
 1      0.81624      0.81624    0.0005%
 2      1.36047      1.36039    0.0055%
 3      1.90522      1.90455    0.0354%
 4      2.45272      2.44871    0.1640%
 5      3.00994      2.99286    0.5706%
```
오차가 1% 미만입니다! 그리고 4.4절에서 얘기했던 것처럼 위쪽 준위로 올라갈수록 오차가 커집니다.  

파동함수도 확인해볼까요?
```python
fig, axes = plt.subplots(1, 2, figsize=(12, 5))
V = 0.5 * k * x**2 / e

axes[0].plot(x * 1e12, V, color="black", lw=2)
axes[1].plot(x * 1e12, V, color="black", lw=2)
for v in range(5):
    f = psi[:, v] / np.sqrt(dx)
    if f[np.argmax(np.abs(f))] < 0:
        f = -f
    axes[0].hlines(E[v], -50, 50, color="gray", lw=1, ls=":")
    axes[1].hlines(E[v], -50, 50, color="gray", lw=1, ls=":")
    axes[0].plot(x * 1e12, E[v] + 0.2 * f / np.abs(f).max(), lw=2)
    axes[1].plot(x * 1e12, E[v] + 0.2 * f**2 / np.abs(f**2).max(), lw=2)
axes[0].set_xlim(-40, 40)
axes[0].set_ylim(0, E[5])
axes[0].set_xlabel("Displacement (pm)")
axes[0].set_ylabel("E (eV)")
axes[0].set_title("psi_v")

axes[1].set_xlim(-40, 40)
axes[1].set_ylim(0, E[5])
axes[1].set_xlabel("Displacement (pm)")
axes[1].set_title("|psi_v|^2")

plt.tight_layout()
plt.show()
```
![양자 조화 진동자의 파동함수](/assets/image-62.png)

각 준위의 에너지 높이에 파동함수를 얹어서 그렸습니다. 상자 속 입자에서와 마찬가지로 에너지 준위가 올라가면서 봉우리가 하나씩 늘어나고 마디도 하나씩 늘어나고 있습니다.  
눈에 띄는 점이 두 가지 있네요. 첫 번째는 파동함수가 퍼텐셜 곡선 바깥으로 새어 나가고 있는 것입니다. 그리고 두 번째는 양 끝의 봉우리만 진폭이 조금 크다는 것입니다. 고전 조화 진동자를 생각해보면 자연스러운 현상이기는 합니다. 양 끝에서는 물체의 속력이 느려지니 발견될 확률도 높아지죠. 상자 속 입자에서는 이런 현상이 없었는데 말이죠.


## 바닥은 가우스

바닥상태인 $v=0$을 보면 Hermite 다항식이 1이라서 파동함수가 순수한 Gauss 함수 형태로 나타납니다.

$$
\psi_0(x) = \left(\frac{\alpha}{\pi}\right)^{1/4}e^{-\frac{\alpha x^2}{2}}
$$
 
확인해봅시다.
 
```python
alpha = mu * omega / hbar

g0 = psi[:, 0] / np.sqrt(dx)
if g0[len(g0) // 2] < 0:
    g0 = -g0

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(x * 1e12, g0 / g0.max(), color="crimson", lw=3, label="numerical")
ax.plot(
    x * 1e12,
    np.exp(-alpha * x**2 / 2),
    color="black",
    lw=2,
    ls="--",
    label="Gaussian",
)
ax.set_xlim(-40, 40)
ax.set_xlabel("Displacement (pm)")
ax.set_yticks([])
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()

print(
    f"최대 차이: {np.max(np.abs(g0 / g0.max() - np.exp(-alpha * x**2 / 2))) * 100:.4f}%"
)
```
![가우스 바닥상태](/assets/image-63.png)
```
최대 차이: 0.0270%
```
5.4절에서 가우스 함수 형태의 파동함수를 가지는 상태는 위치-운동량 불확정성이 $\hbar/2$로 최소가 된다고 했었죠. 계산해봅시다.
 
```python
for v in range(4):
    p = psi[:, v] / np.sqrt(dx)
    x_mean = np.trapezoid(p * x * p, x)
    x2_mean = np.trapezoid(p * x**2 * p, x)
    p2_mean = np.trapezoid(p * (-(hbar**2) * (D2(N, dx) @ p)), x)
    sx, sp = np.sqrt(x2_mean - x_mean**2), np.sqrt(p2_mean)
    print(f"v={v}: sigma_x sigma_p = {sx*sp/hbar:.5f} hbar")
```
```
v=0: sigma_x sigma_p = 0.50000 hbar
v=1: sigma_x sigma_p = 1.50001 hbar
v=2: sigma_x sigma_p = 2.50014 hbar
v=3: sigma_x sigma_p = 3.50121 hbar
```
 
바닥상태의 불확정성이 정확히 $0.5\hbar$입니다. 그리고 위쪽 준위는 $(v + \frac{1}{2})\hbar$로 깔끔하게 늘어납니다.(수치 오차는 무시하고요.) 에너지가 $(v + \frac{1}{2})\hbar\omega$였던 것과 같은 구조죠. 자연에 최소 불확정성 상태가 실제로 존재한다는 것을 확인한 셈입니다.
 
 
## 다음 이야기
 
이제 우리는 양자역학적으로 분자의 진동 운동을 설명할 수 있습니다. 이 결과를 실제 분자에 적용할 차례입니다. 다음 절에서는 적외선 분광학에서 이 개념을 어떻게 써먹는지 알아보겠습니다.
 
 
## 확인 문제
 
1. 힘상수 $k$를 4배로 키우면 준위 간격이 몇 배가 될까요? 영점 에너지는요?
2. D₂ 분자(중수소)의 환산질량은 H₂의 2배입니다. 영점 에너지가 어떻게 달라질까요? 코드로 확인해보세요.
3. $v=10$인 상태의 확률밀도를 그려보세요. 위 코드에서 $v$값만 바꿔서 그리면 이상하게 나올 겁니다. $L$값도 같이 바꿔줘야 하는데 왜 그래야 하는지 퍼텐셜 곡선을 보고 생각해보세요.