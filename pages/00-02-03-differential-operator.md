# 0.2.3. 미분 연산자와 행렬

드디어 선형대수와 미적분학이 합쳐지는 0장의 하이라이트입니다. 0.2.의 도입부에서 선형 연산자는 행렬로 쓸 수 있으니 미분 연산자도 행렬로 쓸 수 있다고 했었죠. 그런데 행렬의 성분을 정하려면 먼저 기저를 정해줘야 합니다. 미분의 기저? *생전 처음 보는 단어의 등장*입니다...  

$$
f(x) \longrightarrow
\begin{bmatrix}
f(x_1) \\
f(x_2) \\
\vdots \\
f(x_N)
\end{bmatrix}
\qquad \frac{d}{dx} \longrightarrow
\begin{bmatrix}
D_{11} & D_{12} & \cdots & D_{1N} \\
D_{21} & D_{22} & \cdots & D_{2N} \\
\vdots & \vdots & \ddots & \vdots \\
D_{N1} & D_{N2} & \cdots & D_{NN}
\end{bmatrix}
\qquad f'(x) = \frac{d}{dx}f(x) \longrightarrow
\begin{bmatrix}
f'(x_1) & f'(x_2) & \dots & f'(x_N)
\end{bmatrix}
$$

0.1.5.에서 함수를 격자점들로 쪼갰던 것을 기억하시나요? 이렇게 구간을 $N$개의 점으로 나누면, 함수는 $N$차원 벡터가 되고, 연산자는 $N \times N$ 크기의 행렬이 됩니다. 이것이 바로 기저를 정하는 일이었습니다. 하지만 미분을 그대로 행렬로 만드려면 점을 무한히 많이 찍어야 하고 그러면 행렬의 크기도 무한히 커집니다. 이렇게는 안 되겠군요. 다른 방법을 찾아야겠습니다.  


## 차분의 행렬 표현

그렇다면 미분 대신 차분을 사용하는 건 어떨까요? 0.2.1.에서 컴퓨터에게 미분을 시킬 때 차분을 사용했었습니다. 한번 중심차분을 가지고 행렬을 만들어볼까요? 중심차분 식을 다시 한번 들여다봅시다.

$$
f'(x_i) \approx \frac{f(x_{i+1}) - f(x_{i-1})}{2h}
$$

$i$번째 점에서의 미분값이 양옆의 두 점의 함숫값만 가지고 계산됩니다. 여기서 행렬과 벡터의 곱셈을 떠올려봅시다. 결과 벡터의 $i$번째 성분은 행렬의 $i$번째 행과 입력 벡터의 내적이었습니다. 중심차분의 식을 이런 방식으로 쓴다면 $i$번째 행의 $i-1$번 자리에는 -1을, $i+1$번 자리에는 1을 넣고 나머지는 0으로 채운 뒤 $2h$로 나누면 됩니다.

$$
D = \frac{1}{2h} \begin{bmatrix}
0 & 1 & 0 & 0 & \cdots \\
-1 & 0 & 1 & 0 & \cdots \\
0 & -1 & 0 & 1 & \cdots \\
0 & 0 & -1 & 0 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \ddots \\
\end{bmatrix}
$$

한번 실제로 만들어보겠습니다. `np.eye(N, k=1)`은 대각선에서 한 칸 위에 1을 넣어주고, `np.eye(N, k=-1)`은 대각선에서 한 칸 아래에 1을 넣어줍니다.
```python
import numpy as np

np.set_printoptions(linewidth=120, suppress=True)


def derivative_matrix(N, h):
    return (np.eye(N, k=1) - np.eye(N, k=-1)) / (2 * h)


print(derivative_matrix(6, 1.0))
```
```
[[ 0.   0.5  0.   0.   0.   0. ]
 [-0.5  0.   0.5  0.   0.   0. ]
 [ 0.  -0.5  0.   0.5  0.   0. ]
 [ 0.   0.  -0.5  0.   0.5  0. ]
 [ 0.   0.   0.  -0.5  0.   0.5]
 [ 0.   0.   0.   0.  -0.5  0. ]]
```
이제 이렇게 만든 행렬이 정말로 작동하는지 확인해보겠습니다. $f(x)=\sin(\pi x/L)$이라는 함수에 적용해보겠습니다. 미분하면 $f'(x)=(\pi/L)\cos(\pi x/L)$이 되죠.
```python
import matplotlib.pyplot as plt

L = 1.0
N = 200

x = np.linspace(0, L, N + 2)[1:-1]
h = L / (N + 1)

f = np.sin(np.pi * x / L)

numerical = derivative_matrix(N, h) @ f
analytic = (np.pi / L) * np.cos(np.pi * x / L)
max_error = np.max(np.abs(numerical - analytic))

plt.plot(x, analytic)
plt.plot(x, numerical, "o", markevery=10)
plt.legend(["analytic", "numerical"])
plt.show()

print(f"최대 오차 = {max_error:.3e}")
```
![1차 차분 행렬](/assets/image-15.png)
```
최대 오차 = 1.279e-04
```
실제로 미분한 함수와 거의 똑같은 값이 나오네요! 행렬을 곱했을 뿐인데 미분(엄밀하게는 차분이지만요.)을 할 수 있게 되었습니다.  
여기서 `[1:-1]`로 양 끝점을 뺀 이유는 잠시 뒤에 설명하겠습니다.  

[[TIP]]
만들어진 행렬을 잘 보면 반대칭행렬입니다. 그런데 이 행렬에 허수단위 $i$를 곱하면 에르미트 행렬이 됩니다! 양자역학에서 운동량을 표현하는 연산자($\hat{p}=-i\hbar \frac{d}{dx}$)에는 허수 단위가 붙어 있는데, 반대칭 행렬에 $i$를 곱하면 에르미트 행렬이 되죠. 그래서 운동량 연산자는 에르미트 연산자가 되고, 운동량은 항상 실수로 나타납니다. 자세한 이야기는 6장에서 계속하겠습니다.
[[/TIP]]


## 2차 차분을 행렬로

슈뢰딩거 방정식에는 2차 미분이 들어있습니다. 그럼 2차 차분은 어떻게 행렬로 만들 수 있을까요? 0.2.1.에서 기억하라고 했던 3개의 숫자, 기억나시나요?

$$
f''(x_i) \approx \frac{f(x_{i+1}) - 2f(x_i) + f(x_{i-1})}{h^2}
$$

바로 2차 차분 식에서 함숫값의 계수로 나타났던 $1, -2, 1$이죠! 아까 1차 차분의 행렬을 만든 것과 동일한 방법으로 이 계수들을 행렬에 심으면 됩니다.

$$
D_2 = \frac{1}{h^2} \begin{bmatrix}
-1 & 1 & 0 & 0 & \cdots \\
1 & -2 & 1 & 0 & \cdots \\
0 & 1 & -2 & 1 & \cdots \\
0 & 0 & 1 & -2 & \cdots \\
\vdots & \vdots & \vdots & \vdots & \ddots \\
\end{bmatrix}
$$

대각선과 그 한 칸 위, 아래에만 값이 있고 나머지는 전부 0입니다. 이렇게 대각선 모양으로만 0이 아닌 값이 채워진 행렬을 **띠행렬(band matrix)**이라 부르고, 그 중에서도 $D_2$ 같이 3줄만 채워진 것을 **삼중대각 행렬(tridiagonal matrix)**이라고 부릅니다.  
```python
def second_derivative_matrix(N, h):
    return (np.eye(N, k=1) - 2 * np.eye(N) + np.eye(N, k=-1)) / (h**2)


print(second_derivative_matrix(6, 1.0))
```
```
[[-2.  1.  0.  0.  0.  0.]
 [ 1. -2.  1.  0.  0.  0.]
 [ 0.  1. -2.  1.  0.  0.]
 [ 0.  0.  1. -2.  1.  0.]
 [ 0.  0.  0.  1. -2.  1.]
 [ 0.  0.  0.  0.  1. -2.]]
```

차분 행렬처럼 대부분의 값이 0인 행렬을 **희소행렬(sparse matrix)**이라고 부릅니다. 저장 공간과 계싼 시간을 크게 줄일 수 있어서 실제 계산에서 매우 유리합니다.
