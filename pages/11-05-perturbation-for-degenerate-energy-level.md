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

