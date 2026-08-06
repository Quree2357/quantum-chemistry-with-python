# 3.3. Schrödinger 방정식

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/03-03.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

드디어 전자의 파동 방정식이 어떤 형태인지 알아볼 때가 되었습니다. 바로 Schrödinger 방정식이죠.  

앞의 두 절에서 우리는 고전적인 파동방정식을 1차원과 2차원에서 풀어봤습니다. 그러면 전자에게도 그 식을 쓰면 되지 않을까요?  
결론부터 말하자면, 안 됩니다. 그리고 왜 안 되는지를 따져보는 과정에서, Schrödinger 방정식이 왜 그런 모양이어야 하는지가 드러납니다.


## 물질파와 파동방정식

2장에서 배운 것을 정리해봅시다. 전자는 파장 $\lambda = h/p$인 파동이고, 에너지는 $E = h\nu$였습니다. 이걸 다시 쓰면 이렇게 됩니다.

$$
p = \hbar k, \qquad E = \hbar\omega
$$

여기서 $k = 2\pi/\lambda$이고 파수(wavenumber)라고 부르며, $\omega = 2\pi\nu$는 각진동수(angular frequency)라고 부릅니다. $2\pi$를 묶어두면 식이 깔끔해져서 물리학에서는 이 형태를 즐겨 씁니다.

이제 가장 단순한 파동 하나를 생각해봅시다. 앞에서 봤던 사인/코사인 함수죠. 지수함수꼴로 쓰면 이렇게 됩니다.

$$
\Psi(x, t) = e^{i(kx - \omega t)}
$$

0.1.1에서 봤던 그 형태죠. 이걸 앞에서 봤던 고전 파동방정식에 넣으면 어떻게 될까요?

$$
\frac{\partial^2 \Psi}{\partial x^2} = \frac{1}{v^2}\frac{\partial^2 \Psi}{\partial t^2}
$$

각각 미분해보면 왼쪽은 $-k^2\Psi$, 오른쪽은 $-\omega^2\Psi/v^2$이 됩니다. 양변을 $\Psi$로 나누면

$$
k^2 = \frac{\omega^2}{v^2} \quad \Longrightarrow \quad \omega = vk
$$

$\omega$가 $k$에 비례합니다. 이런 관계를 **분산 관계(dispersion relation)**라고 부르는데, 방정식이 어떤 파동을 허용하는지 알려주는 지문 같은 것입니다.

그런데 전자는 어떤 분산 관계를 가져야 할까요? 자유 전자의 에너지는 운동에너지뿐이고 우리가 잘 알고 있는 운동에너지 식 $E=\frac{1}{2}mv^2$를 변형하면 다음을 얻습니다.

$$
E = \frac{p^2}{2m}
$$

여기에 $E = \hbar\omega$와 $p = \hbar k$를 넣으면

$$
\hbar\omega = \frac{\hbar^2 k^2}{2m} \quad \Longrightarrow \quad \omega = \frac{\hbar k^2}{2m}
$$

$\omega$가 $k$의 제곱에 비례합니다. 두 관계가 다르죠.

$$
\text{고전 파동: } \omega \propto k \qquad \qquad \text{물질파: } \omega \propto k^2
$$

고전 파동방정식은 전자의 파동을 만들어낼 수 없습니다. 틀린 파동을 내놓거든요.


## Schrödinger 방정식

그러면 $\omega \propto k^2$이 나오는 방정식을 한번 만들어봅시다.

미분이 무슨 일을 하는지 보세요. $\Psi = e^{i(kx-\omega t)}$에 대해 $x$로 한 번 미분하면 $ik$가 딸려 나오고, $t$로 한 번 미분하면 $-i\omega$가 딸려 나옵니다. 그러니 $\omega$의 1차와 $k$의 2차를 얻으려면 시간에 대해서는 한 번, 공간에 대해서는 두 번 미분해야 합니다.

$$
\frac{\partial \Psi}{\partial t} = C\frac{\partial^2 \Psi}{\partial x^2}
$$

$C$가 어떤 상수인지 알아봅시다. $\Psi$를 대입해보면 다음과 같죠.

$$
-i\omega \Psi = C(-k^2 \Psi) \quad \Longrightarrow \quad \omega = -iCk^2
$$

여기서 중요한 일이 벌어집니다. $\omega$는 진동수이니 실수여야 합니다. 그런데 우변에 $-i$가 붙어 있죠. 그러니 $C$ 자체가 순허수여야 $\omega$가 실수로 나옵니다.

허수 단위가 방정식 안으로 들어올 수밖에 없는 겁니다.

$C = i\hbar/2m$으로 두면 $\omega = \hbar k^2/2m$이 정확히 나옵니다. 정리하면 이렇게 되겠네요.

$$
i\hbar\frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \Psi}{\partial x^2}
$$

확인해봅시다. SymPy로 자유 전자 파동을 넣어보죠.

```python
import sympy as sp

x, t = sp.symbols("x t")
k, w, m, hbar = sp.symbols("k omega m hbar", positive=True)

psi = sp.exp(sp.I * (k * x - w * t))

left = sp.I * hbar * sp.diff(psi, t)
right = -(hbar**2) / (2 * m) * sp.diff(psi, x, 2)

print("좌변 / psi =", sp.simplify(left / psi))
print("우변 / psi =", sp.simplify(right / psi))
```
```
좌변 / psi = hbar*omega
우변 / psi = hbar**2*k**2/(2*m)
```

좌변이 $\hbar\omega = E$이고 우변이 $\frac{\hbar^2k^2}{2m} = \frac{p^2}{2m}$입니다. 방정식이 말하는 바가 곧 $E = \frac{p^2}{2m}$이었던 것이죠.

만약 전자가 퍼텐셜 장 안에서 움직인다면 퍼텐셜 에너지까지 더해주면 됩니다.

$$
i\hbar\frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \Psi}{\partial x^2} + V(x)\Psi
$$

이것이 바로 **시간 의존 Schrödinger 방정식(time-dependent Schrödinger equation)**입니다. 1926년에 발표되었죠.

읽는 법은 간단합니다. 우변의 첫 항이 운동에너지, 둘째 항이 퍼텐셜 에너지에 해당하는 항입니다.  결국 에너지 보존 법칙을 파동의 언어로 쓴 것입니다.


## 그래서 무엇이 다른가

고전 파동방정식과 나란히 놓아봅시다.

$$
\frac{\partial^2 u}{\partial x^2} = \frac{1}{v^2}\frac{\partial^2 u}{\partial t^2}
\qquad\qquad
i\hbar\frac{\partial \Psi}{\partial t} = -\frac{\hbar^2}{2m}\frac{\partial^2 \Psi}{\partial x^2} + V \Psi
$$

- 첫째, 시간 미분이 1차입니다.  
고전 파동방정식은 시간에 대해 2차죠. 이 차이 때문에 초기조건도 달라집니다. 고전 파동방정식에서는 처음 모양과 처음 속도를 둘 다 알아야 하지만, 전자는 처음 파동함수 하나만 주어지면 이후가 전부 정해집니다.

- 둘째, 허수 단위 $i$가 식 안에 있습니다.  
이건 편의를 위해 붙인 게 아닙니다. 방금 본 것처럼 $i$가 없으면 $\omega$가 실수로 나오지 않죠.

- 셋째, 해가 일반적으로 복소수입니다.  
고전 파동에서는 복소수가 계산 요령이었습니다. 0.1.1에서 "실수부만 꺼내면 원래 코사인으로 돌아온다"고 했었죠. 그런데 여기서는 실수부만 꺼내면 방정식을 만족하지 않습니다.

[[TIP]]
$i$가 왜 필요한지, 0.2.3에서 이미 복선이 있었습니다. 그때 1차 미분 행렬 $D_1$이 반대칭이라 `eigh`를 쓸 수 없다고 했죠. 그리고 반대칭 행렬에 $i$를 곱하면 에르미트 행렬이 된다고요.
지금 상황이 정확히 그것입니다. 시간 미분이 1차라서 반대칭적인 성격을 갖는데, 여기에 $i$를 곱해야 에르미트가 되고 그래야 에너지가 실수로 나옵니다. 6장에서 이 이야기를 다시 제대로 하겠습니다.
[[/TIP]]


## 이 차이가 만드는 결과

분산 관계가 다르면 실제로 무슨 일이 벌어질까요? 2.4절에서 만든 파동 묶음을 시간에 따라 움직여봅시다.

```python
import numpy as np
import matplotlib.pyplot as plt

x = np.linspace(-40, 90, 3000)
k0, sk = 2.0, 0.25
k = np.linspace(k0 - 6 * sk, k0 + 6 * sk, 1200)
A = np.exp(-((k - k0) ** 2) / (2 * sk**2))


def packet(t, dispersion):
    w = dispersion(k)
    return (A[:, None] * np.exp(1j * (k[:, None] * x[None, :] - w[:, None] * t))).sum(
        axis=0
    )


classical = lambda kk: 1.0 * kk  # omega = v k
quantum = lambda kk: 0.5 * kk**2  # omega = hbar k^2 / 2m  (hbar = m = 1)

fig, axes = plt.subplots(2, 3, figsize=(15, 6), sharey="row")

for j, t in enumerate([0, 15, 30]):
    for i, (disp, name, col) in enumerate(
        [
            (classical, "classical wave", "steelblue"),
            (quantum, "matter wave", "crimson"),
        ]
    ):
        p = np.abs(packet(t, disp)) ** 2
        axes[i, j].plot(x, p / p.max(), color=col, lw=1.8)
        axes[i, j].set_xlim(-20, 90)
        axes[i, j].set_yticks([])
        axes[i, j].set_title(f"{name},  t = {t}", fontsize=10)
        axes[i, j].set_xlabel("x")

plt.tight_layout()
plt.show()
```
![슈뢰딩거 방정식](/assets/image-38.png)

위쪽 파란 곡선이 고전 파동입니다. 모양을 그대로 유지하면서 오른쪽으로 움직입니다. 3.1절에서 $f(x-vt)$가 해라고 했던 것과 같죠. 반면에 아래쪽 빨간 곡선은 물질파인데, 움직이면서 점점 퍼지는 것을 볼 수 있습니다.  
폭을 재보면 이렇습니다.

```python
def width(t, disp):
    p = np.abs(packet(t, disp)) ** 2
    p = p / np.trapezoid(p, x)
    mu = np.trapezoid(x * p, x)
    return np.sqrt(np.trapezoid((x - mu) ** 2 * p, x))


print(f"{'t':>4}  {'고전파 폭':>10}  {'물질파 폭':>10}")
for t in [0, 10, 20, 30]:
    print(f"{t:4d}  {width(t, classical):10.3f}  {width(t, quantum):10.3f}")
```
```
   t    고전파 폭    물질파 폭
   0       2.828       2.828
  10       2.828       3.335
  20       2.828       4.528
  30       2.828       6.010
```

고전파는 폭이 그대로인데 물질파는 두 배 이상 퍼졌습니다.

이유는 분산 관계에 있습니다. 고전 파동에서는 $\omega = vk$이니 파장이 달라도 속도가 전부 같습니다. 그래서 여러 파장이 섞여 있어도 다 같이 움직여 모양이 유지되죠. 반면 물질파는 $\omega \propto k^2$이라 파장마다 속도가 다릅니다. 빠른 성분이 앞서가고 느린 성분이 처지면서 파동 패킷이 벌어지는 겁니다.

그리고 이건 2.4절의 이야기와 이어집니다. 파동 묶음이 퍼진다는 것은 위치의 불확정성이 커진다는 뜻이니까요. 전자를 한 곳에 가둬놓아도 시간이 지나면 어디 있는지 점점 모르게 됩니다.


## Schrödinger 방정식은 양자역학의 공리다

여기서 분명히 해둘 것이 있습니다. 방금 우리가 한 것은 Schrödinger 방정식을 유도한 게 아닙니다.

"분산 관계가 이래야 하니 미분 횟수가 이렇고, 그러니 $i$가 필요하다"는 이야기는 그럴듯하지만, 어디까지나 왜 이런 모양인지 납득해보는 과정입니다. 실제로 우리는 $E=p^2/2m$이라는 답을 이미 알고 거꾸로 맞춰 들어갔죠. McQuarrie와 같은 다른 양자역학 책들은 다른 출발점에서 시작해서 Schrödinger 방정식을 이끌어내기도 하지만 마찬가지로 유도한 것은 아닙니다. Schrödinger 방정식은 더 근본적인 무언가에서 유도된 것이 아닙니다. Newton의 $F=ma$가 유도된 식이 아닌 것과 같습니다. 근거는 실험 결과와 잘 맞는다는 것 하나뿐이죠. 그래서 과학자들은 Schrödinger 방정식을 양자역학의 공리 중 하나로서 받아들였습니다.

그리고 100년 가까이 지난 지금까지, 이 방정식이 틀렸다는 실험 결과는 하나도 나오지 않았습니다. 마치 Newton의 운동법칙이 틀렸다는 실험 결과가 하나도 없는 것처럼요.

[[TIP]]
사실 Schrödinger 방정식과 Newton의 운동법칙 모두 상대론적인 효과를 무시한 상태에서 세워진 이론이라 속도가 매우 빨라 광속에 가까워지는 경우에서는 정확해지지 않습니다. 하지만 지금 여기서는 크게 신경쓰지 않아도 됩니다.
[[/TIP]]


## 시간을 독립시키기

Schrödinger 방정식은 시간 변수를 포함하고 있는데, 사실 우리가 앞으로 다룰 문제 대부분은 시간에 따라 변하지 않는 상황입니다. 그러면 시간 부분을 따로 떼어내기 위해 3.1절에서 썼던 변수분리법을 그대로 적용해봅시다. 상수는 $E$라고 놓읍시다. 왜냐면 나중에 이 값이 에너지 자리에 놓일 거거든요. 만약 퍼텐셜 $V$가 시간에 무관하다면 $\Psi(x,t) = \psi(x)T(t)$로 둘 수 있고, 정리하면 시간 부분이 이렇게 떨어져 나옵니다.

$$
T(t) = e^{-iEt/\hbar}
$$

그리고 공간 부분은 시간 변수가 사라진 방정식을 만족하게 됩니다.

$$
-\frac{\hbar^2}{2m}\frac{d^2\psi}{dx^2} + V(x)\psi = E\psi
$$

이것을 **시간 독립 Schrödinger 방정식(time-independent Schrödinger equation)**이라고 부릅니다. 앞으로 이 책에서 "Schrödinger 방정식"이라고 하면 대부분 이것을 뜻합니다. 그리고 앞으로 대문자 $\Psi$가 아닌 소문자 $\psi$라고 하면 시간 독립인 경우의 해를 의미합니다.

여기서 재미있는 점이 있습니다. 시간 부분 $e^{-iEt/\hbar}$는 크기가 1인 복소수라서, 확률밀도를 구하면 통째로 사라집니다.

$$
|\Psi(x,t)|^2 = |\psi(x)|^2 \cdot |e^{-iEt/\hbar}|^2 = |\psi(x)|^2
$$

확률밀도가 시간에 따라 변하지 않습니다. 파동함수 자체는 계속 위상이 돌고 있는데, 우리가 관측하는 양은 멈춰 있는 것이죠. 이런 상태를 **정상상태(stationary state)**라고 부릅니다.

1.3절이 떠오르시나요? Bohr가 근거 없이 가정했던 것 중 하나가 "특정 궤도에 있는 동안에는 빛을 내지 않는다"였습니다. 그 궤도를 정상상태라고 불렀고요. 그 정상상태의 수학적 정체가 바로 이것입니다. 확률 분포가 시간에 따라 변하지 않으니 고전적인 의미에서의 전자기파 방출은 나타나지 않습니다. 2.2절에서는 각운동량에 대한 가정을 갚았고, 여기서는 두 번째 가정을 갚은 셈입니다.

이제 우리는 Schrödinger 방정식을 손에 넣었습니다. 그런데 이 방정식을 만족한다고 아무 함수나 전자의 파동함수가 될 수 있는 건 아닙니다. 다음 절에서는 파동함수가 갖춰야 할 조건들을 정리하겠습니다.


## 확인 문제

1. 자유 전자의 분산 관계 $\omega = \hbar k^2/2m$에서 위상속도 $v_p = \omega/k$를 구해보세요. 이 값이 전자의 실제 속도 $p/m$과 같나요?
2. 위 시뮬레이션에서 `sk`(파수 폭)를 절반으로 줄여보세요. 파동 묶음이 퍼지는 속도가 어떻게 달라지나요? 2.4절의 불확정성 원리로 설명할 수 있을까요?
3. 고전 파동방정식의 해 $\cos(kx-\omega t)$는 실수입니다. 이 함수를 Schrödinger 방정식에 넣으면 어떻게 될까요? 왜 해가 될 수 없는지 설명해보세요.
