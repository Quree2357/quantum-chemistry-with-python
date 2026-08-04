# 0.2.2. 수치적분

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/00-02-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

0.1.5에서 함수의 내적을 계산할 때 `np.trapezoid`라는 함수를 아무 설명 없이 썼습니다. 이 절에서 그 정체를 밝히겠습니다.

이미 적분과 Riemann 합 이야기를 했으니 사실 절반은 아시는 셈입니다. 다만 "잘게 나눠서 더한다"는 아이디어에도 여러 가지 방법이 있습니다.


## 직사각형으로 근사하기

가장 단순한 방법부터 봅시다. 구간을 잘게 나누고, 각 조각을 직사각형으로 채우는 겁니다. 높이는 왼쪽 끝 함숫값으로 잡고요.

$$
\int_a^b f(x)\,dx \approx \sum_i f(x_i)\,h
$$

0.1.5에서 `np.dot(f, g) * dx`로 계산했던 게 바로 이것입니다. Riemann 합의 정의 그대로죠.

앞 절과 마찬가지로 답을 아는 함수로 시험해봅시다. $f(x) = xe^x$를 $[0,1]$에서 적분하면 정확히 1이 나옵니다. 0.1.5에서 SymPy로 확인했던 그 함수입니다.

```python
import numpy as np

f = lambda x: x * np.exp(x)
exact = 1.0

x = np.linspace(0, 1, 5)
h = x[1] - x[0]
y = f(x)

rect = np.sum(y[:-1]) * h
print(f"직사각형 : {rect:.10f}   오차 = {abs(rect - exact) / abs(exact) *100:.2f} %")
```
```
직사각형 : 0.6832792505   오차 = 31.67 %
```

오차가 32%나 됩니다... 왜 그런지는 그림으로 보면 바로 알 수 있습니다.

```python
import matplotlib.pyplot as plt

xf = np.linspace(0, 1, 300)
xn = np.linspace(0, 1, 5)

fig, ax = plt.subplots(figsize=(5.5, 4))
for i in range(len(xn) - 1):
    ax.add_patch(
        plt.Rectangle(
            (xn[i], 0),
            xn[i + 1] - xn[i],
            f(xn[i]),
            fc="steelblue",
            ec="navy",
            alpha=0.45,
        )
    )
ax.plot(xf, f(xf), color="black", lw=2, zorder=5)
ax.plot(xn, f(xn), "o", color="black", ms=5, zorder=6)
ax.set_xlabel("x")
ax.grid(alpha=0.3)
plt.show()
```
![직사각형 근사](/assets/image-11.png)

함수는 계속 올라가는데 계산하는 넓이는 직사각형으로 잘라버렸습니다. 그래서 매 조각마다 곡선과 직사각형 사이에 빈틈이 생기죠. 그 빈틈들이 전부 오차입니다.


## 사다리꼴로 근사하기

그러면 직사각형 말고 사다리꼴을 쓰면 어떨까요? 조각의 왼쪽 끝과 오른쪽 끝을 직선으로 이어주는 겁니다. 위쪽을 비스듬히 잘라내니 훨씬 잘 맞겠죠.

$$
\int_a^b f(x)\,dx \approx \sum_i \frac{h}{2}\left(f(x_i)+f(x_{i+1})\right)
$$

이게 **사다리꼴 공식(trapezoidal rule)**이고, `np.trapezoid`가 하는 계산입니다.

```python
trap = np.trapezoid(y, x)
print(f"사다리꼴 : {trap:.10f}   오차 = {abs(trap - exact) / abs(exact) *100:.2f} %")
```
```
사다리꼴 : 1.0230644791   오차 = 2.31 %
```

똑같이 5개의 점을 썼는데 오차가 거의 10배가 줄었습니다. 그림으로 다시 봅시다.

```python
fig, ax = plt.subplots(figsize=(5.5, 4))
for i in range(len(xn) - 1):
    ax.fill(
        [xn[i], xn[i + 1], xn[i + 1], xn[i]],
        [0, 0, f(xn[i + 1]), f(xn[i])],
        fc="crimson",
        ec="darkred",
        alpha=0.4,
    )
ax.plot(xf, f(xf), color="black", lw=2, zorder=5)
ax.plot(xn, f(xn), "o", color="black", ms=5, zorder=6)
ax.set_xlabel("x")
ax.grid(alpha=0.3)
plt.show()
```
![사다리꼴 근사](/assets/image-12.png)

아까 있던 삼각형 빈틈이 거의 사라졌습니다. 사다리꼴의 윗변이 실제 곡선보다 위에 있어 여전히 아주 살짝 어긋나긴 하지만, 직사각형과는 비교가 안 되죠.


## 포물선으로 근사하기

여기서 한 발 더 나가봅시다. 사다리꼴은 두 점을 직선으로 이었습니다. 직선 대신 곡선으로 이으면 더 좋지 않을까요? Simpson 공식은 세 점을 포물선으로 이어서 넓이를 구합니다. SciPy에서 사용해볼 수 있습니다.

$$
\int_a^b f(x)\,dx \approx \sum_i \frac{x_{i+1}-x_i}{6}\left[f(x_i) + 4f\left(\frac{x_i + x_{i+1}}{2}\right) + f(x_{i+1})\right]
$$

```python
from scipy.integrate import simpson

simp = simpson(y, x=x)
print(f"Simpson   : {simp:.10f}   오차 = {abs(simp - exact) / abs(exact) * 100:.3f} %")
```
```
Simpson   : 1.0001690471   오차 = 0.017 %
```

또 한 번 크게 좋아졌습니다. 사다리꼴보다 약 100배쯤 더 정확하네요.

```python
fig, ax = plt.subplots(figsize=(5.5, 4))
for i in range(0, len(xn) - 1, 2):
    xs = np.array([xn[i], xn[i + 1], xn[i + 2]])
    c = np.polyfit(xs, f(xs), 2)
    xx = np.linspace(xs[0], xs[-1], 80)
    ax.fill_between(
        xx, 0, np.polyval(c, xx), fc="darkorange", ec="chocolate", alpha=0.4
    )
    ax.plot(xx, np.polyval(c, xx), color="chocolate", lw=1.6, ls="--", zorder=4)
ax.plot(xf, f(xf), color="black", lw=2, zorder=5)
ax.plot(xn, f(xn), "o", color="black", ms=5, zorder=6)
ax.set_xlabel("x")
ax.grid(alpha=0.3)
plt.show()
```
![Simpson 근사](/assets/image-13.png)

주황 점선이 세 점을 지나는 포물선입니다. 검은 곡선과 거의 겹쳐서 구분이 안 될 정도네요. 점 5개만으로 이 정도라니 놀랍습니다.

[[TIP]]
예전 SciPy에서는 `scipy.integrate.simps`라는 이름이었지만 `simpson`으로 바뀌었습니다. `np.trapz`가 `np.trapezoid`로 바뀐 것과 비슷한 상황이니, 오래된 코드를 참고하실 때 주의하세요.
[[/TIP]]


## 세 방법을 나란히 놓고 보면
 
이제 셋을 한눈에 비교해봅시다. 위에서 그린 그림들을 나란히 놓고, 각 방법의 오차도 함께 표시하겠습니다.
 
```python
y = f(xn)
h = xn[1] - xn[0]
errs = [
    abs(np.sum(y[:-1]) * h - exact),
    abs(np.trapezoid(y, xn) - exact),
    abs(simpson(y, x=xn) - exact),
]
titles = ["rectangle", "trapezoid", "Simpson"]

fig, axes = plt.subplots(1, 3, figsize=(12, 3.8), sharey=True)

for k, (ax, title) in enumerate(zip(axes, titles)):
    if k == 0:
        for i in range(len(xn) - 1):
            ax.add_patch(
                plt.Rectangle(
                    (xn[i], 0),
                    xn[i + 1] - xn[i],
                    f(xn[i]),
                    fc="steelblue",
                    ec="navy",
                    alpha=0.45,
                )
            )
    elif k == 1:
        for i in range(len(xn) - 1):
            ax.fill(
                [xn[i], xn[i + 1], xn[i + 1], xn[i]],
                [0, 0, f(xn[i + 1]), f(xn[i])],
                fc="crimson",
                ec="darkred",
                alpha=0.4,
            )
    else:
        for i in range(0, len(xn) - 1, 2):
            xs = np.array([xn[i], xn[i + 1], xn[i + 2]])
            c = np.polyfit(xs, f(xs), 2)
            xx = np.linspace(xs[0], xs[-1], 80)
            ax.fill_between(
                xx, 0, np.polyval(c, xx), fc="darkorange", ec="chocolate", alpha=0.4
            )
    ax.plot(xf, f(xf), color="black", lw=2, zorder=5)
    ax.plot(xn, f(xn), "o", color="black", ms=5, zorder=6)
    ax.set_xlabel("x")
    ax.grid(alpha=0.3)
    ax.set_title(f"{title}\nerror = {errs[k]:.2e}", fontsize=11)

plt.show()
```
![근사법 비교](/assets/image-14.png)
 
계단처럼 뭉텅뭉텅 잘라내는 직사각형, 비스듬히 따라가는 사다리꼴, 곡선에 거의 달라붙는 Simpson. 같은 점 5개를 쓰는데도 오차는 $10^{-1}$, $10^{-2}$, $10^{-4}$로 자릿수가 통째로 달라집니다.


## 얼마나 빨리 좋아질까요

점의 개수를 늘려가며 세 방법을 비교해봅시다.

```python
for N in [11, 101, 1001, 10001]:
    x = np.linspace(0, 1, N)
    h = x[1] - x[0]
    y = f(x)
    r = abs(np.sum(y[:-1]) * h - exact)
    t = abs(np.trapezoid(y, x) - exact)
    s = abs(simpson(y, x=x) - exact)
    print(f"N ={N:6d}   직사각형 {r:.3e}   사다리꼴 {t:.3e}   Simpson {s:.3e}")
```
```
N =    11   직사각형 1.322e-01   사다리꼴 3.696e-03   Simpson 4.366e-06
N =   101   직사각형 1.355e-02   사다리꼴 3.697e-05   Simpson 4.374e-10
N =  1001   직사각형 1.359e-03   사다리꼴 3.697e-07   Simpson 4.374e-14
N = 10001   직사각형 1.359e-04   사다리꼴 3.697e-09   Simpson 0.000e+00
```

점을 10배로 늘렸을 때 오차가 줄어드는 비율을 보세요.

| 방법 | 오차 감소 | 정확도 |
|---|---|---|
| 직사각형 | 10분의 1 | $O(h)$ |
| 사다리꼴 | 100분의 1 | $O(h^2)$ |
| Simpson | 10000분의 1 | $O(h^4)$ |

앞 절에서 전진차분이 $O(h)$, 중심차분이 $O(h^2)$였던 것과 같은 이야기입니다. 근사를 조금만 영리하게 하면 정확도의 **차수 자체**가 올라가죠.

Simpson은 $N=10001$에서는 오차가 아예 0으로 찍혔습니다. 물론 진짜 0이라서 그런 건 아니고, 컴퓨터가 표현할 수 있는 한계보다 오차가 작아진 것입니다.


## 적분은 잘게 나눌수록 좋아집니다

앞 절에서 미분은 $h$를 너무 줄이면 오히려 나빠졌습니다. 적분도 그럴까요?

```python
for N in [10**3, 10**4, 10**5, 10**6]:
    x = np.linspace(0, 1, N)
    err = abs(np.trapezoid(f(x), x) - exact)
    print(f"N = {N:9d}   오차 = {err:.3e}")
```
```
N =      1000   오차 = 3.705e-07
N =     10000   오차 = 3.698e-09
N =    100000   오차 = 3.697e-11
N =   1000000   오차 = 3.697e-13
```

정확도가 계속 좋아집니다. 왜 그럴까요? 앞 절에서 미분이 나빠졌던 이유는 거의 같은 두 수를 빼면서 앞자리가 상쇄되어 사라졌기 때문이었습니다. 수치적분도 반올림 오차에서 완전히 자유롭지는 않지만, 미분처럼 거의 같은 두 값을 빼고 작은 \(h\)로 나누는 과정이 없기 때문에, 일반적으로는 격자를 상당히 촘촘하게 만들 때까지 정확도가 안정적으로 좋아집니다.

그래서 수치계산에서는 이런 말이 있습니다. *적분은 쉽고 미분은 어렵다.* 수학 수업에서 배운 것과 정반대라 재미있죠. 손으로 할 때는 미분이 기계적이고 적분이 까다로웠는데, 컴퓨터에서는 뒤집힙니다.

물론 공짜는 아닙니다. $N=10^7$이면 계산 시간이 그만큼 오래 걸리니까요. 정확도와 비용을 저울질하는 이 문제는 앞으로 계속 만나게 됩니다.


## 그래서 이걸 어디에 써요?

앞으로 이 책에서 적분은 지겹도록 나옵니다. 대표적인 두 가지만 미리 보여드리죠.

첫째, **정규화(normalization)**입니다. 파동함수는 $\int|\psi|^2dx = 1$을 만족해야 한다고 3장에서 배우게 됩니다. 확률의 총합이 1이라는 뜻이죠. 0.1.5에서 만들었던 사인 함수로 확인해봅시다.

```python
L = 1.0
x = np.linspace(0, L, 2001)
psi = np.sin(np.pi * x / L)

norm2 = np.trapezoid(psi**2, x)
A = 1 / np.sqrt(norm2)

print(f"적분값       = {norm2}")
print(f"규격화 상수 A = {A}")
print(f"sqrt(2/L)   = {np.sqrt(2/L)}")
```
```
적분값       = 0.5
규격화 상수 A = 1.414213562373095
sqrt(2/L)   = 1.4142135623730951
```

0.1.5에서 기저함수 앞에 붙였던 $\sqrt{2/L}$의 정체가 이것이었습니다. 그때는 그냥 크기를 1로 맞추는 계수라고만 했었죠.

둘째, **기댓값(expectation value)**입니다. 6장에서 어떤 물리량의 평균값을 구할 때 $\int\psi^*\hat{A}\psi\,dx$ 같은 적분을 계산하게 됩니다. 전자의 평균 위치, 평균 운동에너지 같은 것들이죠.

그리고 19장에 가면 이런 적분을 수십만 번 계산하게 됩니다. 그때가 되면 적분 하나의 비용이 왜 중요한지 뼈저리게 느끼시게 될 겁니다.


## 정리하며

이제 미분과 적분을 모두 컴퓨터로 계산할 수 있게 되었습니다. 다음 절에서는 이 절의 진짜 목적지로 갑니다. 미분을 행렬로 만들 차례입니다.


## 확인 문제

1. $\int_0^\pi \sin x\,dx = 2$입니다. 세 방법으로 각각 계산해 오차를 비교해보세요.
2. Simpson 공식은 왜 짝수 개의 구간(홀수 개의 점)을 요구할까요? 세 점씩 묶어 포물선을 그린다는 사실로부터 생각해보세요.
3. 0.1.5에서 계산했던 $\langle e_1, e_1 \rangle = 1$을 `simpson`으로 다시 계산해보세요. `trapezoid`보다 정확한가요?
