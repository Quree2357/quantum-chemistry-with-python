# 11.1. 변분 원리

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/11-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

먼저 **변분법(variational method)**에 대해서 살펴봅시다. Euler와 Lagrange에 의해 정립된 변분법은 양자역학뿐만 아니라 수학과 물리학 전반에서 쓰이는 방법인데, 함수의 최대/최소값을 찾는 방식입니다. 5.3절에서 이런 이야기를 잠깐 하고 넘어갔었죠.

> 아무 함수나 가져와 에너지 기댓값을 계산해도 이 부등식이 성립합니다. (...) 이 성질을 변분 원리라고 부르는데, (...) 11장에서 제대로 다루겠습니다.

변분법은 정답을 몰라도 근사의 정확도를 채점할 수 있다는 점에서 매우 유용한 방법입니다. 한번 자세히 들여다보죠.


## 정답을 모르면 상한이라도 알자

9장에서 본 것처럼 헬륨 원자부터는 Schrödinger 방정식을 정확하게 풀 수 있는 방법이 없습니다. 생각해볼 수 있는 방법 중 하나는 여러 모양의 파동 함수를 추측해서 방정식에 던져보고 해가 되는지 확인하는 것인데, 문제는 그게 정말 정답인지 알 수가 없다는 것입니다. 정답이 없으니 비교해 볼 수도 없고 정답이 있었다면 애초에 추측할 필요가 없으니까요. 이럴 때 **변분 원리(variational principle)**가 이 악의 무리의 순환 공격을 끊어주는 역할을 해줍니다.  
그냥 '변분법'에 대해서 검색해보면 미적분학과 고전역학에서 사용하는 내용들이 나올 텐데 지금 우리가 사용할 변분 원리는 양자역학에서 사용할 형태입니다. 임의의 계(system)의 바닥 상태를 생각해봅시다. 바닥 상태에서의 파동 함수가 $\psi_0$이고 그때의 에너지가 $E_0$라고 해보죠. 그러면 변분 원리는 우리에게 아무거나 함수를 하나 골라서 에너지를 계산해보면 항상 $E_0$보다 크거나 같은 값을 준다는 것을 알려줍니다. 그러니까 추측을 여러 개 던져봤을 때 그 중에서 제일 낮은 에너지를 가지는 것이 정답에 가깝다는 뜻입니다. 정답에 얼마나 가까운지는 알 수 없어도 적어도 어떤 경우가 제일 정확도가 높은지는 알 수 있는 거죠. 식으로 쓰면 이렇게 됩니다.

$$
E_{\phi} = \frac{\left\langle \phi \middle| \hat{H} \middle| \phi \right\rangle}{\langle \phi | \phi \rangle} \geq E_0
$$

여기서 우리가 아무렇게나 잡은 함수 $\phi$를 **시험 함수(trial function)**라고 부릅니다. 증명 과정은 아주 간단합니다. 5.3절에서 절반 정도는 해놨는데요, 시험 함수 $\phi$를 고유 함수들의 선형 결합으로 나타낸 것에서 시작합니다.

$$
\phi = \sum_n c_n \psi_n
$$

이 때, $\psi_n$은 Hamiltonian의 정확한 고유 함수들입니다. 즉, $\hat{H}\psi_n = E_n \psi_n$을 만족하는 함수들이죠. 이것을 위의 에너지 기댓값 식에 대입하면 다음과 같습니다. (시험 함수 $\phi$가 규격화되어 있다고 가정합시다.)

$$
E_{\phi} = \left\langle \phi \middle| \hat{H} \middle| \phi \right\rangle = \int{\phi^* \hat{H} \phi} d\tau = \int{ \sum_n c_n^* \psi_n^* \hat{H} \sum_n c_n \psi_n } d\tau = \int{ \sum_n c_n^* \psi_n^* \sum_n c_n \hat{H} \psi_n } d\tau
$$

여기서 $\psi_n$은 정확한 고유 함수들이니 고유값 방정식을 만족했었죠. 그러면 이렇게 됩니다.

$$
E_{\phi} = \int{ \left( \sum_n c_n^* \psi_n^* \right) \left( \sum_m c_m E_m \psi_m \right) } d\tau = \sum_n \sum_m c_n^* c_m E_m \int{\psi_n^* \psi_m} \, d\tau
$$

고유함수의 직교성을 이용하면 위 식에서 $n=m$일 때만 남고 나머지는 모두 0이 되겠죠.

$$
E_{\phi} = \sum_n c_n^* c_n E_n = \sum_n |c_n|^2 E_n
$$

이제 5.3절에서 봤던 값이 나왔습니다. 바닥 상태의 에너지가 $E_n$ 중에서 제일 낮은 값일 테니 $E_{\phi}$는 항상 $E_0$보다 크거나 같습니다. 변분 원리가 증명되었습니다!

[[TIP]]
등호는 당연히 시험 함수가 정확한 바닥 상태의 파동 함수일 때 성립합니다. 그러니 에너지를 낮추려는 방향으로 계산하면 자연스럽게 바닥 상태에 가까워지는 거죠. 다만 주의할 점이 있습니다. 에너지가 잘 맞는다고 파동 함수도 잘 맞는 것은 아닙니다. 에너지가 소수점 아래 몇 자리까지 맞아도 쌍극자 모멘트 같은 다른 성질은 몇 퍼센트씩 틀리는 일이 생깁니다. 이 이야기는 책의 후반부에서 다시 하겠습니다.
[[/TIP]]


## 답을 아는 문제로 검증

우리가 지금까지 풀었던 문제 중에 답을 정확히 아는 문제에 변분 원리를 적용해보겠습니다. 상자 속 입자 문제가 적당해보이네요. 4.4절에서 계산했던 바로는 $L=1$ nm 상자 속 입자의 바닥 상태 에너지는 0.376 eV였습니다.  
시험 함수는 경계 조건만 만족하면 정말로 아무거나 넣어봐도 됩니다. 상자 속 입자 문제에서는 경계 조건이 양 끝에서 0이 되는 것이었죠. 가장 단순한 함수는 포물선입니다.

$$
\phi = x(L-x) \qquad (0 \leq x \leq L)
$$

상자 안에서는 $V=0$이니 Hamiltonian은 $T$ 하나만으로 구성되죠.

$$
\left\langle \phi \middle| \hat{H} \middle| \phi \right\rangle = \left\langle \phi \middle| \hat{T} \middle| \phi \right\rangle = -\frac{\hbar^2}{2m} \int_0^L{\phi^* \, \frac{d^2}{dx^2} \phi} \, dx
$$

우리가 정한 함수는 규격화되어 있지 않으니 분모도 계산해주어야 합니다. SymPy로 직접 계산해봅시다.
```python
import sympy as sp
from scipy.constants import hbar, m_e, e, h

L = 1e-9
E_exact = h**2 / (8 * m_e * L**2) / e
print(f"참값 : {E_exact:.6f} eV")

x, L_, hbar_, m = sp.symbols("x L hbar m", positive=True)

psi1 = x * (L_ - x)
psi2 = x**2 * (L_ - x)**2
psi3 = sp.sin(sp.pi * x / L_)

num1 = sp.integrate(-hbar_**2 / (2 * m) * psi1 * psi1.diff(x, 2), (x, 0, L_))
denum1 = sp.integrate(psi1**2, (x, 0, L_))
energy1 = num1 / denum1
energy1_eV = energy1.evalf(subs={hbar_: hbar, L_: L, m: m_e}) / e
print(f"x(L-x): {energy1} = {energy1_eV:.6f} eV (오차: {(energy1_eV - E_exact) / E_exact * 100:+8.4f} %)")

num2 = sp.integrate(-hbar_**2 / (2 * m) * psi2 * psi2.diff(x, 2), (x, 0, L_))
denum2 = sp.integrate(psi2**2, (x, 0, L_))
energy2 = num2 / denum2
energy2_eV = energy2.evalf(subs={hbar_: hbar, L_: L, m: m_e}) / e
print(f"x²(L-x)²: {energy2} = {energy2_eV:.6f} eV (오차: {(energy2_eV - E_exact) / E_exact * 100:+8.4f} %)")

num3 = sp.integrate(-hbar_**2 / (2 * m) * psi3 * psi3.diff(x, 2), (x, 0, L_))
denum3 = sp.integrate(psi3**2, (x, 0, L_))
energy3 = num3 / denum3
energy3_eV = energy3.evalf(subs={hbar_: hbar, L_: L, m: m_e}) / e
print(f"sin(πx/L): {energy3} = {energy3_eV:.6f} eV (오차: {(energy3_eV - E_exact) / E_exact * 100:+8.4f} %)")
```
```
참값 : 0.376030 eV
x(L-x): 5*hbar**2/(L**2*m) = 0.380998 eV (오차:  +1.3212 %)
x²(L-x)²: 6*hbar**2/(L**2*m) = 0.457198 eV (오차: +21.5854 %)
sin(πx/L): pi**2*hbar**2/(2*L**2*m) = 0.376030 eV (오차:  +0.0000 %)
```

정확한 해석해인 사인 함수와 비교해봤습니다. 확실히 변분 원리대로 시험 함수는 참값보다 큰 에너지를 주네요. 그런데 조금 더 복잡한 함수인 $x^2(L-x)^2$은 오차가 매우 크게 납니다. 왜 그럴까요? 함수의 모양을 보면 알 수 있습니다.
```python
import numpy as np
import matplotlib.pyplot as plt

xs = np.linspace(0, L, 500)
fig, ax = plt.subplots(figsize=(7, 4))
for f, name, col in [
    (xs * (L - xs), "x(L-x)", "crimson"),
    (xs**2 * (L - xs) ** 2, "x²(L-x)²", "steelblue"),
    (np.sin(np.pi * xs / L), "sin(πx/L)", "black"),
]:
    ax.plot(xs * 1e9, f / f.max(), lw=2, color=col, label=name)
ax.set_xlabel("x (nm)")
ax.set_ylabel("psi (normalized)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![시험 함수의 비교](/assets/image-86.png)

4차 포물선 함수는 2차 포물선 함수에 비해서 사인 함수와 너무 안 맞습니다. 특히 양쪽 끝에서 차이가 크게 나죠. 일반적으로 다항식의 차수가 높아질수록 근사가 잘 되는 건 맞지만 지금처럼 모양이 안 맞는 함수를 넣으면 오히려 오차가 더 커진다는 것을 알 수 있습니다. 추가로 4차 포물선 함수로 근사한 경우도 비교해보죠. 이번에는 적분이 복잡하게 나올테니 결과값만 출력하도록 하겠습니다.
```python
psi4 = x * (L_ - x) * sp.pi / L_**2 + x**2 * (L_ - x)**2 * (16 - 4 * sp.pi) / L_**4

num4 = sp.integrate(-hbar_**2 / (2 * m) * psi4 * psi4.diff(x, 2), (x, 0, L_))
denum4 = sp.integrate(psi4**2, (x, 0, L_))
energy4 = num4 / denum4
energy4_eV = energy4.evalf(subs={hbar_: hbar, L_: L, m: m_e}) / e
print(f"4차 근사: {energy4_eV:.6f} eV (오차: {(energy4_eV - E_exact) / E_exact * 100:+8.4f} %)")
```
```
4차 근사: 0.376040 eV (오차:  +0.0026 %)
```
오차가 매우 작아져서 소수점 넷째 자리까지 맞네요. 이런 식으로 정확한 파동 함수의 모양과 점점 비슷한 시험 함수를 넣으면 참값에 가까워지는 값을 얻게 됩니다.


## 매개변수와 변분법

그렇다고 정답에 가까운 값이 나올 거라고 믿으면서 무턱대고 아무 시험 함수나 막 집어넣어볼 수는 없습니다. 어떻게 하면 효율적으로 시험 함수를 찾을 수 있을까요? 한 가지 방법은 시험 함수에 매개변수를 넣어놓고 값을 조절하면서 에너지가 가장 낮아지도록 만드는 값을 찾는 것입니다. 이번에는 조화 진동자 문제에 이 방법을 적용해보죠.  
조화 진동자의 퍼텐셜을 보면 원점에서는 0이고 거리가 멀어질수록 급격하게 무한대로 발산합니다. 그러니 파동 함수의 경계 조건은 거리가 무한대일 때 0이 되는 것이죠. 그리고 원점 기준으로 대칭이어야 할테니 $e^{-x^2}$ 같은 함수가 적당해 보이네요. 길이의 -2승 차원을 갖는 매개변수 $c>0$를 지수에 넣어서 시험 함수 $\phi = e^{-cx^2}$를 만들겠습니다.

$$
\left\langle \hat{H} \right\rangle = -\frac{\hbar^2}{2m} \int_{-\infty}^{\infty}{\phi^* \, \frac{d^2}{dx^2} \phi} \, dx + \frac{1}{2}m\omega^2 \int_{-\infty}^{\infty}{\phi^* x^2 \phi } \, dx \qquad \phi = e^{-cx^2}
$$

이 시험 함수는 규격화되지 않은 상태이기 때문에 $\langle \phi | \phi \rangle$로 나누어주어야 합니다. SymPy로 계산해봅시다.
```python
from sympy import oo
c, w = sp.symbols("c omega", positive=True)
psi5 = sp.exp(-c * x**2)

num5 = sp.integrate(-hbar_**2 / (2 * m) * psi5 * psi5.diff(x, 2) + m * w**2 * psi5 * x**2 * psi5 / 2, (x, -oo, oo))
denum5 = sp.integrate(psi5**2, (x, -oo, oo))
print(sp.simplify(num5 / denum5))
```
```
c*hbar**2/(2*m) + m*omega**2/(8*c)
```

$$
E_{\phi} = \frac{\hbar^2 c}{2m} + \frac{m \omega^2}{8c}
$$

이제 에너지가 매개변수 $c$에 대한 식으로 나타났습니다. 에너지의 최소값을 찾는 방법은 단순하게 이 식을 $c$에 대해 미분해서 0이 되는 값을 찾아주면 되죠.

$$
\frac{dE}{dc} = \frac{\hbar^2}{2m} - \frac{m \omega^2}{8c^2}=0
$$

이 식을 풀면 우리가 찾는 $c$의 값은 $\frac{m \omega}{2\hbar}$임을 알 수 있습니다. 이 값을 다시 에너지 식에 대입해보면 $E=\frac{1}{2}\hbar\omega$가 나옵니다. 우리가 6.2절에서 구했던 값과 정확히 일치하는 값이죠. 이것은 우리가 시험 함수의 형태를 정확한 바닥 상태 파동 함수와 똑같이 잡았기 때문입니다.

[[TIP]]
들뜬 상태에 대해서도 변분법을 쓸 수 있을까요? 사실 위에서 살펴본 변분법은 바닥 상태에 대해서만 적용 가능합니다. 다만 조건을 하나 붙이면 되는데요. 시험 함수가 바닥 상태의 파동 함수와 직교하도록 정하면 됩니다. 그러면 $c_0$이 0이 되면서 $E_0$ 항이 빠지기 때문에 그 다음 준위에 대해서 같은 논리가 성립하거든요. 하지만 실제로는 바닥 상태의 파동 함수가 정확히 어떤 형태인지 모르니 이 조건을 그대로 가져가기 어렵습니다. 그래서 들뜬 상태에 대한 계산이 바닥 상태보다 훨씬 까다롭죠.
[[/TIP]]


## 다음 이야기

매개변수를 포함한 시험 함수의 모양을 정하고 에너지를 계산한 후에 에너지가 최소가 되는 매개변수의 값을 찾는 것, 이것이 변분법의 전부입니다. 다음 절에서는 변분법을 조금 더 깊게 들어가서 매개변수를 여러 개 사용하는 경우에 대해 살펴보겠습니다.


## 확인 문제
1. 상자 속 입자 코드에서 시험 함수를 $\sin^2(\pi x/L)$로 넣어서 에너지를 계산해보세요. 오차가 어느 정도 나나요?
2. 상자 속 입자 문제에서 시험 함수의 매개변수를 $[x(L-x)]^c$처럼 넣어서 최적화해보세요. $c$의 값이 얼마일 때 에너지가 최소가 되나요?
3. 어떤 계산 결과가 실험값보다 낮게 나왔다면 무엇을 의심해야 할까요? 변분 원리의 관점에서 생각해보세요.
