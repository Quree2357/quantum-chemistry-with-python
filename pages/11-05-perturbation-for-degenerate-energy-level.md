# 11.5. 축퇴 준위의 섭동

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/11-05.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

11.3절에서 섭동론은 축퇴가 없는 경우에만 계산이 가능했습니다. 축퇴가 있으면 에너지가 같은 상태가 존재하니 보정 식에서 분모가 0이 되어 발산해버리거든요. 그런데 축퇴는 생각보다 흔하게 일어납니다. 4장에서 봤던 3차원 상자 속 입자나 8장의 수소 원자 등의 경우에 모두 축퇴가 존재했죠. 이번 절에서는 축퇴가 있을 때 섭동론을 어떻게 적용할 수 있는지 살펴보겠습니다.


## 축퇴가 있으면 상태가 섞인다

만약에 어떤 계의 에너지 준위가 $d$-중으로 축퇴되어 있다고 가정합시다. 이 때의 에너지를 $E^{(0)}$라고 하고 축퇴되어 있는 상태들의 파동 함수를 $\psi_i^{(0)}$이라고 해보죠. 그러니까 다음과 같은 상황입니다.

$$
\hat{H}^{(0)}\psi_i^{(0)} = E^{(0)}\psi_i^{(0)} \qquad (i=1,2,\dots,d)
$$

이건 지금 섭동이 없는 상태의 Schrödinger 방정식이죠. 이제 여기에 섭동을 가해야 합니다.  

$$
\hat{H} = \hat{H}^{(0)} + \lambda \hat{H}'
$$

그런데 다시 말하지만 축퇴가 있으면 보정 항을 계산할 때 분모가 0이 되는 문제가 발생합니다. 이걸 어떻게 해결하면 될까요?  

잠깐 생각해보면, 보정 항에서 분자는 다른 상태들이 얼마나 섞이는지를 나타낸다고 했습니다. 그러니까 분모가 0에 가까워진다는 것은 상태들이 엄청나게 크게 섞인다는 뜻이죠. 여기서 아이디어를 얻어봅시다. 아예 처음부터 상태들을 섞어놓고 출발점으로 삼아보자는 것입니다. 이렇게요.

$$
\phi_n^{(0)} = \sum_{i=1}^d c_i \psi_i^{(0)} \qquad (n=1,2,\dots,d)
$$

각 $n$에 대해서 계수들의 값이 달라질 겁니다. 이제 이 계수들이 얼마인지 알면 됩니다.

다음 단계는 비축퇴 섭동론에서 본 것과 비슷한 과정입니다. 에너지와 파동 함수를 $\lambda$에 대해 전개하는 거죠. 대신에 $\psi_n^{(0)}$ 자리에 $\phi_n^{(0)}$를 넣습니다.

$$
\psi_n = \phi_n^{(0)} + \lambda \psi_n^{(1)} + \lambda^2 \psi_n^{(2)} + \cdots
$$
$$
E_n = E^{(0)} + \lambda E_n^{(1)} + \lambda^2 E_n^{(2)} + \cdots
$$

이제 이걸 원래 방정식에 대입합니다. 

$$
\left( \hat{H}^{(0)} + \lambda \hat{H}' \right) \left( \phi_n^{(0)} + \lambda \psi_n^{(1)} + \cdots \right) = \left( E^{(0)} + \lambda E_n^{(1)} + \cdots \right) \left( \phi_n^{(0)} + \lambda \psi_n^{(1)} + \cdots \right)
$$

이걸 모두 다 전개할 필요는 없고요, $\lambda$의 차수가 같은 항들의 계수 비교만 하면 됩니다. 먼저 0차부터 볼까요?

$$
\hat{H}^{(0)} \phi_n^{(0)} = E^{(0)} \phi_n^{(0)}
$$

아... 시시해서 그만 죽...을 정도까지는 아니고 원래 만족해야 하는 식이 나왔으니 그냥 패스하면 됩니다. 다음 1차 항을 봅시다.

$$
\left( \hat{H}^{(0)} - E^{(0)} \right) \psi_n^{(1)} = \left( E_n^{(1)} - \hat{H}' \right) \phi_n^{(0)}
$$

정리하면 이렇게 나오는데요. 이제 양변에 $\psi_m^{(0)}$를 곱하고 적분해봅시다.

$$
\left\langle \psi_m^{(0)} \middle| \hat{H}^{(0)} - E^{(0)} \middle| \psi_n^{(1)} \right\rangle = E_n^{(1)} \left\langle \psi_m^{(0)} \middle| \phi_n^{(0)} \right\rangle - \left\langle \psi_m^{(0)} \middle| \hat{H}' \middle| \phi_n^{(0)} \right\rangle
$$

여기에서 주목할 점은 $\hat{H}^{(0)}$ 연산자는 Hermitian이고, $\psi_m^{(0)}$은 고유값 $E^{(0)}$에 대응하는 고유 함수이니 첫 번째 적분이 0이 된다는 것입니다.

$$
\left\langle \psi_m^{(0)} \middle| \hat{H}' \middle| \phi_n^{(0)} \right\rangle - E_n^{(1)} \left\langle \psi_m^{(0)} \middle| \phi_n^{(0)} \right\rangle = 0
$$

이제 처음에 선형 결합으로 썼던 $\phi_n^{(0)}$를 대입해봅시다.

$$
\sum_{i=1}^d c_i \left\langle \psi_m^{(0)} \middle| \hat{H}' \middle| \psi_i^{(0)} \right\rangle - E_n^{(1)} \sum_{i=1}^d c_i \left\langle \psi_m^{(0)} \middle| \psi_i^{(0)} \right\rangle = 0
$$

우리는 여기서 '사실 $\psi_m^{(0)}$들은 원래 모두 정규직교한다'라고 할 수 있습니다! 이게 하이라이트입니다. 그러면 두 번째 항의 적분이 Kronecker delta가 됩니다.

$$
\sum_{i=1}^d \left[ \left\langle \psi_m^{(0)} \middle| \hat{H}' \middle| \psi_i^{(0)} \right\rangle - E_n^{(1)} \delta_{mi} \right] c_i = 0 \qquad m=1,2,\dots,d
$$

이제 $d$개의 식을 가진 연립방정식을 풀기만 하면 됩니다. 어? 그러면 행렬로 쓸 수 있겠는데요?

$$
\begin{bmatrix}
\hat{H}'_{11} & \hat{H}'_{12} & \cdots & \hat{H}'_{1d} \\
\hat{H}'_{21} & \hat{H}'_{22} & \cdots & \hat{H}'_{2d} \\
\vdots & \vdots & \ddots & \vdots \\
\hat{H}'_{d1} & \hat{H}'_{d2} & \cdots & \hat{H}'_{dd} 
\end{bmatrix} \mathbf{c} = 
\begin{bmatrix}
E_n^{(1)} \\
E_n^{(2)} \\
\vdots \\
E_n^{(d)}
\end{bmatrix} \mathbf{c}
$$

어떤 교재는 왼쪽에 있는 행렬을 **섭동 행렬(perturbation matrix)**이라고 부르고 $\mathbf{W}$라고 쓰기도 하는데요, 이 쪽이 쓰기 편하니 가져다 쓰겠습니다. 우리는 이 방정식이 해를 가질 조건을 알고 있습니다.

$$
\det (\mathbf{W} - \mathbf{E} \cdot \mathbf{I}) = 0
$$

이 방정식을 풀면 근이 $d$개가 나올 겁니다. 각각이 비섭동 축퇴 에너지 준위에 대한 1차 보정 값이죠. 만약 모든 근이 서로 다른 값을 가진다면 축퇴가 완전히 깨져서 $d$개의 상태로 갈라질 거고요. 그런데 가끔 중근이 나오는 경우도 있습니다. 그러면 축퇴가 완전히 깨지지는 않고 더 적은 수의 상태로 갈라질 겁니다.  

와, 정말 길었습니다. 이제는 실제로 적용해볼 시간입니다.


## 찌그러진 상자 속 입자

4.3절로 돌아가봅시다. 이번에는 상자를 약간 찌그러뜨려서 정육면체에서 살짝 달라진 경우를 살펴보죠. $b=(1+\varepsilon)a$와 $c=(1+2\varepsilon)a$로 정해볼까요? 그러면 에너지가 이렇게 됩니다.

$$
E = E_0 \left( n_x^2 + \frac{n_y^2}{(1+\varepsilon)^2} + \frac{n_z^2}{(1+2\varepsilon)^2} \right)
$$

$\varepsilon$이 아주 작다고 가정하면 $\frac{1}{(1+\varepsilon)^2} \approx 1-2\varepsilon$으로 근사할 수 있으니 다음과 같이 쓸 수 있습니다.

$$
E \approx E_0 (n_x^2 + n_y^2 + n_z^2) - \varepsilon E_0 (2n_y^2 + 4n_z^2)
$$

뒤의 항이 섭동이 되겠네요. 이 경우에는 정말 운이 좋은 것이 섭동 행렬이 대각 행렬이 됩니다. 서로 다른 상태들이 섞이지 않는 거죠.

```python
import numpy as np

eps = 0.02
states = [(1, 1, 2), (1, 2, 1), (2, 1, 1)]

W = np.diag([-eps * (2 * ny**2 + 4 * nz**2) for (nx, ny, nz) in states])

print("섭동 행렬 W (E_0 단위):")
print(np.round(W, 4))

E1, C = np.linalg.eigh(W)
print(f"\n1차 보정: {np.round(E1, 4)}")
print(f"보정 에너지: {np.round(6 + E1, 4)}")


def exact(a, b, c, state):
    nx, ny, nz = state
    return (nx / a) ** 2 + (ny / b) ** 2 + (nz / c) ** 2


a, b, c = 1, 1 + eps, 1 + 2 * eps
print(f"\n{'상태':>8} {'1차 보정':>9} {'정확한 값':>8}")
for s in states:
    approx = 6 - eps * (2 * s[1] ** 2 + 4 * s[2] ** 2)
    print(f"{str(s):>10} {approx:10.4f} {exact(a, b, c, s):11.4f}")
```
```
섭동 행렬 W (E_0 단위):
[[-0.36  0.    0.  ]
 [ 0.   -0.24  0.  ]
 [ 0.    0.   -0.12]]

1차 보정: [-0.36 -0.24 -0.12]
보정 에너지: [5.64 5.76 5.88]

      상태     1차 보정    정확한 값
 (1, 1, 2)     5.6400      5.6594
 (1, 2, 1)     5.7600      5.7692
 (2, 1, 1)     5.8800      5.8857
```
세 상태가 갈라졌고 순서도 일치합니다. 섭동이 작아서 1차 보정만으로도 에너지 값이 어느 정도 맞네요.

그런데 이 예시에서는 섭동이 상태를 섞지 않았습니다. 축퇴 섭동론의 진짜 힘을 아직 보여주지 못한 것이죠.


## 수소 원자의 Stark 효과

수소 원자에 자기장을 걸면 Zeeman 효과에 의해 축퇴가 깨진다고 했었는데요, 비슷한 일이 전기장을 걸어도 나타납니다. 이 현상을 발견한 독일의 물리학자의 이름을 따서 Stark 효과라고 부릅니다. (프리렌에 나오는 슈타르크가 아닙니다!)

만약 전기장이 $z$ 방향으로 걸려 있다고 해보죠. 그러면 전기장에 의해 전자에 추가적으로 퍼텐셜 에너지 $\Phi(z) = eEz$가 걸립니다. 이 부분이 섭동이죠. 전기장이 충분히 약하다면 섭동론을 적용할 수 있습니다.

수소 원자의 $n=1$ 준위에는 축퇴가 없으니 축퇴가 있는 $n=2$ 준위를 살펴보죠. $s$ 오비탈과 $p$ 오비탈 3개로 총 4중 축퇴가 되어 있는 상태입니다. 섭동 행렬을 만들어봅시다.

먼저 대각 성분부터 볼까요?

$$
\left\langle \psi_{2s} \middle| eEz \middle| \psi_{2s} \right\rangle, \quad \left\langle \psi_{2p} \middle| eEz \middle| \psi_{2p} \right\rangle
$$

오비탈의 모양을 생각해보면 파동 함수가 $xy$-평면에 대해 대칭이니 여기에 기함수인 $z$를 곱해서 적분하면 0이 된다는 사실을 금방 알 수 있습니다.

대각 성분이 아닌 성분들 중에서는 $\left\langle \psi_{2s} \middle| eEz \middle| \psi_{2p_z} \right\rangle$만 0이 아닙니다. SymPy로 확인해봅시다.

```python
import sympy as sp

r, theta, phi = sp.symbols("r theta phi", positive=True)

psi_2s = 1 / (4 * sp.sqrt(2 * sp.pi)) * (2 - r) * sp.exp(-r / 2)
psi_2pz = 1 / (4 * sp.sqrt(2 * sp.pi)) * r * sp.cos(theta) * sp.exp(-r / 2)
z = r * sp.cos(theta)


def braket(f, op, g):
    integrand = f * op * g * r**2 * sp.sin(theta)
    return sp.integrate(sp.integrate(sp.integrate(integrand, (phi, 0, 2 * sp.pi)), (theta, 0, sp.pi)), (r, 0, sp.oo))

print("<2s|z|2s>   =", sp.simplify(braket(psi_2s, z, psi_2s)))
print("<2pz|z|2pz> =", sp.simplify(braket(psi_2pz, z, psi_2pz)))
print("<2s|z|2pz>  =", sp.simplify(braket(psi_2s, z, psi_2pz)))
```
```
<2s|z|2s>   = 0
<2pz|z|2pz> = 0
<2s|z|2pz>  = -3
```
$2p_x$와 $2p_y$는 $z$축에 수직한 방향으로 뻗어있어서 $z$와 섞이지 않아 0이 됩니다. 그러니 섭동 행렬은 이렇게 되겠네요.

$$
\mathbf{W} = eE
\begin{bmatrix}
0 & -3 & 0 & 0 \\
-3 & 0 & 0 & 0 \\
0 & 0 & 0 & 0 \\
0 & 0 & 0 & 0
\end{bmatrix}
$$

순서는 $(2s, 2p_z, 2p_x, 2p_y)$로 뒀습니다. 이제 이 행렬을 대각화해봅시다.

```python
eE = 0.01

W = eE * np.array([[0, -3, 0, 0], [-3, 0, 0, 0], [0, 0, 0, 0], [0, 0, 0, 0]], dtype=float)

E1, C = np.linalg.eigh(W)

print("1차 보정:", np.round(E1, 5))
print("고유벡터 (열 방향):")
print(np.round(C, 4))
```
```
1차 보정: [-0.03  0.    0.    0.03]
고유벡터 (열 방향):
[[-0.7071  0.      0.     -0.7071]
 [-0.7071  0.      0.      0.7071]
 [-0.      1.      0.      0.    ]
 [-0.      0.      1.      0.    ]]
```
4중 축퇴였던 것이 3중 축퇴(-0.03, 0, 0.03)로 갈라졌습니다. 고유벡터를 보면 갈라진 두 상태가 $2s$와 $2p_z$를 절반씩 섞은 것입니다.

$$
\psi_{\pm} = \frac{1}{\sqrt{2}} \left( \psi_{2s} \mp \psi_{2p_z} \right)
$$

나머지 $2p_x$와 $2p_y$는 그대로 남죠.

섭동 행렬에 $eE$가 곱해져 있으니 전기장의 세기가 커질수록 이 갈라짐도 커질 겁니다. 특히 이 경우에는 선형으로 비례하기 때문에 **1차 Stark 효과(First-order Stark effect)**라고 부릅니다.


## 상태가 섞인다는 것

그런데 $2s$와 $2p_z$가 섞였다는 것이 무슨 뜻일까요? 전기장이 걸리면 수소 원자의 퍼텐셜은 더 이상 순수한 중심력장이 아닙니다. $z$ 방향이 특별해지니까요. 그러면 $\hat{L}^2$과 Hamiltonian이 교환하지 않고, $l$이 확정된 상태가 고유 상태가 아니게 됩니다. 위에서 봤던 섞인 상태는 $l=0$도 아니고 $l=1$도 아닌 거죠.

물리적으로는 두 오비탈이 섞여서 전자 구름이 한 쪽으로 쏠리게 됩니다. 전기장이 걸려있으니 당연히 그렇게 되겠죠. 원자가 전기장에 반응해 쌍극자를 만드는 것입니다.

[[TIP]]
수소 원자에서 1차 Stark 효과가 나타나는 것은 사실 특이한 일입니다. 대부분의 원자에서는 전기장에 대한 에너지의 변화가 제곱으로 비례하거든요. 그 이유는 우연 축퇴 때문입니다. 수소 원자에서는 $2s$와 $2p$가 같은 에너지인데 반전 대칭성이 반대라서 섞일 수 있습니다. 하지만 다른 원자에서는 $2s$와 $2p$ 오비탈이 다른 에너지를 가지니 축퇴가 없어져서 비축퇴 섭동론을 쓰게 됩니다. 이 경우에는 1차 보정이 0이 되서 2차부터 기여하는데, 이걸 2차 Stark 효과라고 부릅니다.
[[/TIP]]


## 정리하며

이 장에서 근사법 두 가지를 배웠습니다. 변분법과 섭동론이죠. 각각의 방법을 어떻게 적용하는지, 그리고 장단점은 무엇인지 봤습니다. 그리고 축퇴가 있을 때 섭동론에서는 선형 변분법이 쓰이는 것도 있었죠. 앞으로도 이 도구들이 계속 쓰일 겁니다.
다음 장부터는 Part III의 시작입니다. 원자에 대한 것은 얼추 봤으니 이제 분자로 넘어가야겠죠. 원자들이 결합한 것이 분자인데 이제는 원자핵도 개수가 늘어납니다. 더욱 복잡해질 것만 같네요. 잠시 재정비를 하고 상급 몬스터들을 잡으러 떠나봅시다!


## 확인 문제
1. 위 상자 코드에서 `eps`를 0.1로 바꿔보세요. 1차 보정값의 오차가 어떻게 되나요?
2. 수소 원자의 $n=3$ 준위는 9중 축퇴입니다. 전기장을 걸면 몇 개로 갈라질까요?
3. 1차 Stark 효과의 갈라짐 $3eEa_0$을 eV 단위로 바꿔보세요. 만약 전기장의 크기가 $10^5$ V/m 정도면 얼마나 갈라질까요?
