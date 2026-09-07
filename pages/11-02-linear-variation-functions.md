# 11.2. 선형 변분 함수

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/11-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

변분법을 사용하면 정답을 몰라도 에너지가 낮아지는 방향을 찾아가면 그게 정답에 가까워지는 길이라는 것을 알 수 있었습니다. 매개변수를 사용해서 시험 함수를 조절해가며 에너지를 최소화하는 방법도 봤고요. 그런데 매개변수를 여러 개 쓴다면 좀 더 에너지가 최소가 되는 지점을 잘 찾을 수 있을 것 같습니다. 매개변수를 하나의 차원에서만 조절하면 다른 방향은 알 수 없으니까요. 물론 그렇게 하면 여러 개의 연립방정식을 풀어야 하니 좀 더 복잡해질 겁니다. 매개변수를 10개 쓰면 10개의 식을 연립해서 풀어야 하니까요. 하지만 매개변수를 넣는 방식을 잘 고르면 이 문제를 해결할 수 있습니다.


## 변분 함수의 선형 결합

시험 함수로 넣을 수 있는 함수는 경계 조건만 만족하면 아무거나 가능하다고 했었죠. 그러니까 이런 시험 함수들을 한꺼번에 다 넣어볼 수도 있을 겁니다. 특히 함수들을 선형 결합 형태로 합친 것을 많이 쓰고, **선형 변분 함수(linear variation function)**라고 부릅니다. 그리고 각각의 시험 함수 $f_i$들은 **기저 함수(basis function)**라고 부릅니다.

$$
\phi = c_1 f_1 + c_2 f_2 + \dots + c_n f_n = \sum_{i=1}^n c_i f_i
$$

이제 이 함수를 앞에서 봤던 에너지 계산식에 넣어보죠. 먼저 규격화를 위한 분모 부분부터 봅시다.

$$
\langle \phi | \phi \rangle = \int{\phi^* \phi} \, d\tau = \int{\left( \sum_{i=1}^n c_i^* f_i^* \right) \left( \sum_{j=1}^n c_j f_j \right)} \, d\tau = \sum_{i=1}^n \sum_{j=1}^n c_i^* c_j \int{f_i^* f_j} \, d\tau = \sum_{i=1}^n \sum_{j=1}^n c_i^* c_j S_{ij}
$$

$S_{ij}=\int{f_i^* f_j} \, d\tau$는 **중첩 행렬(overlap matrix) 요소**라고 하는 값인데, 어디서 많이 본 것처럼 생겼지만 시험 함수들이 직교한다는 보장이 없기 때문에 Kronecker delta와 같지는 않습니다. 이번에는 분자 부분을 살펴봅시다.

$$
\langle \phi | \hat{H} | \phi \rangle = \int{\phi^* \hat{H} \phi} \, d\tau = \int{\left( \sum_{i=1}^n c_i^* f_i^* \right) \hat{H} \left( \sum_{j=1}^n c_j f_j \right)} \, d\tau = \sum_{i=1}^n \sum_{j=1}^n c_i^* c_j \int{f_i^* \hat{H} f_j} \, d\tau = \sum_{i=1}^n \sum_{j=1}^n c_i^* c_j H_{ij}
$$

Hamiltonian이 붙은 것 빼고는 아까와 비슷한 형태입니다. 이번에도 뒤에 붙은 적분을 $H_{ij}$라고 썼는데 Hamiltonian 행렬 요소라고 부릅니다. 이제 에너지는 다음과 같이 계산할 수 있습니다.

$$
E_{\phi} = \frac{\sum_{ij} c_i^* c_j H_{ij}}{\sum_{ij} c_i^* c_j S_{ij}}
$$

이제 에너지의 최소값을 찾기 위해 이 식을 각 $c_i$에 대해 편미분한 값을 0으로 놓고 연립방정식을 풀면 됩니다...라고 하면 말로는 간단해보이지만 사실 방정식 $n$개를 연립해야 하죠. 다행히 잘 정리하면 $n$개의 식이 모두 일차 방정식으로 나옵니다. 중간 과정은 생략하고 결과만 보면 다음과 같습니다.

$$
\sum_{k=1}^n \left[ \left( H_{ik} - S_{ik} E_{\phi} \right) c_k \right]=0 \qquad i=1,2,\dots,n
$$

이걸 어떻게 풀지 싶으시겠지만 우리에게는 연립 일차 방정식을 푸는 엄청난 도구가 있죠. 바로 행렬입니다! 위 연립방정식을 행렬 형태로 쓰면 다음과 같습니다.

$$
\mathbf{H} \mathbf{c} = E_{\phi} \mathbf{S} \mathbf{c}
$$

아하, 고유값 문제였네요! 오른쪽에 $\mathbf{S}$가 붙어있다는 점만 빼면요. 이런 식을 일반화된 고유값 문제라고 합니다. 이제 우리는 이 문제를 `eigh`를 사용해서 풀 수 있습니다!  
참고로 위 식을 손으로 풀려면 우변을 왼쪽으로 넘겨서 방정식 $\det \left( \mathbf{H} - E_{\phi} \mathbf{S} \right) = 0$를 풀면 되는데, 이 식에는 **영년 방정식(secular equation)**이라는 이름이 붙어있습니다.  

아무튼 우리는 이 문제를 풀면 $n$개의 해(고유값)를 얻게 되는데, 이 중에서 가장 작은 값이 바닥 상태 에너지의 상한이 됩니다. 나머지는 각각의 들뜬 상태의 에너지에 대한 상한이 되고요. 그리고 각 에너지에 대응하는 고유벡터 $\mathbf{c}$를 구하면 그 상태의 파동 함수 $\phi$ 를 구성할 수 있게 됩니다.


## 이번에도 상자 속 입자

마찬가지로 상자 속 입자 문제에 이 방법을 적용해봅시다. 시험 함수들을 $f_k=[x(L-x)]^k$라고 잡아보겠습니다.($k$는 자연수입니다.)
```python
import numpy as np
from scipy.constants import hbar, m_e, e, h
from scipy.linalg import eigh

L = 1e-9
x = np.linspace(0, L, 100001)
E_exact = h**2 / (8 * m_e * L**2) / e


def basis(n):
    basis_function = []
    for k in range(1, n + 1):
        f = (x * (L - x)) ** k
        basis_function.append(f / np.sqrt(np.trapezoid(f**2, x)))
    return basis_function


def matrices(B):
    S = np.array([[np.trapezoid(a * b, x) for b in B] for a in B])
    H = np.array([[hbar**2 / (2 * m_e) * np.trapezoid(np.gradient(a, x) * np.gradient(b, x), x) for b in B] for a in B])
    return H, S


B = basis(3)
H, S = matrices(B)
E, C = eigh(H, S)

print(f"참값 (n = 1): {E_exact:.10f} eV")
print(f"계산 결과    : {E[0] / e:.10f} eV")
print(f"오차        : {(E[0] / e - E_exact) / E_exact * 100:+.8f} % ")
```
```
참값 (n = 1): 0.3760301621 eV
계산 결과    : 0.3760301633 eV
오차        : +0.00000031 %
```
기저 함수 3개만으로도 에너지가 소수점 아래 8번째 자리까지 맞습니다! 오차가 사실상 0이나 다름없는 결과입니다.


## 더 많은 기저 함수!

3개만 써도 충분히 참값에 가까운 값을 내놓긴 하지만, 기저 함수의 개수를 늘려가면서 확인해봅시다.
```python
print(f"{'Basis':>4} {'Energy (eV)':>14} {'Error':>14}")

for n in range(1, 5):
    H, S = matrices(basis(n))
    E, C = eigh(H, S)
    err = (E[0] / e - E_exact) / E_exact * 100
    print(f"{n:5d} {E[0] / e:14.8f} {err:+12.6f} %")
```
```
Basis    Energy (eV)          Error
    1     0.38099821    +1.321184 %
    2     0.37603569    +0.001471 %
    3     0.37603016    +0.000000 %
    4     0.37603016    -0.000000 %
```
기저 함수의 개수가 많아질수록 계산되는 에너지가 점점 감소하면서 참값에 가까워지죠. 변분 원리가 제대로 작용하고 있습니다. 선택지가 넓어지니 정확한 파동 함수에 가까워지기가 쉬워지는 것입니다.

[[TIP]]
$n=3$ 이상에서는 에너지 값이 같게 나옵니다. 그 이상으로는 기저 함수의 개수를 늘려도 정확도가 더 늘어나지 않는 것이죠. 실제 계산에서는 이 최적점을 찾아서 계산량과 정확도의 균형을 맞추도록 합니다. 그리고 특히 기저 함수 집합이 너무 커지면 중복 행렬 $\mathbf{S}$가 수치적으로 특이해져서 계산이 아예 이상하게 되는 경우가 있습니다. 이걸 선형 종속성 문제라고 부릅니다.
[[/TIP]]


## 밑에서 진행되고 있는 일

사실 `eigh(H, S)` 함수가 알아서 풀어주기는 하는데, 실제로 안에서 무슨 일이 벌어지는지 알아두면 좋습니다. 기본적인 발상은 기저 함수들이 직교가 되도록 하는 것입니다. 그러면 $\mathbf{S}=\mathbf{I}$가 되면서 보통의 고유값 문제($\mathbf{H} \mathbf{c} = E \mathbf{c}$)로 바뀌거든요. $\mathbf{S}$는 Hermitian이니 0.1.4에서 본 것처럼 대각화할 수 있습니다.

$$
\mathbf{S} = \mathbf{U} \mathbf{s} \mathbf{U}^{\dagger}
$$

$\mathbf{s}$는 고유값들로 이루어진 대각행렬이죠. 그러면 이런 행렬을 만들어볼 수 있습니다.

$$
\mathbf{X} = \mathbf{U} \mathbf{s}^{-1/2} \mathbf{U}^{\dagger}
$$

이렇게 만들어진 행렬 $\mathbf{X}$는 다음과 같은 성질을 갖게 됩니다.

$$
\mathbf{X}^{\dagger} \mathbf{S} \mathbf{X} = \mathbf{I}
$$

이제 기존의 계수 벡터 $\mathbf{c}$를 $\mathbf{X} \mathbf{c}$로 두면 일반화 고유값 문제가 보통의 고유값 문제로 바뀌게 됩니다.

$$
\mathbf{X}^{\dagger} \mathbf{H} \mathbf{X} \mathbf{c} = E \mathbf{c}
$$

직접 계산해봅시다.
```python
H, S = matrices(basis(3))

s, U = np.linalg.eigh(S)
X = U @ np.diag(s**-0.5) @ U.T

print("X^T S X =")
print(np.round(X.T @ S @ X, 10))

H_prime = X.T @ H @ X
E_prime, C_prime = np.linalg.eigh(H_prime)

print(f"직교화 후 에너지: {E_prime[0] / e:.8f} eV")
```
```
X^T S X =
[[ 1. -0.  0.]
 [-0.  1. -0.]
 [ 0. -0.  1.]]
직교화 후 에너지: 0.37603016 eV
```
$\mathbf{X}^{\dagger} \mathbf{S} \mathbf{X} = \mathbf{I}$이 성립하고, 정확히 같은 답이 나왔습니다. 이 방법을 **대칭 직교화(symmetric orthogonalization)**라고 부릅니다. 18장에서 Hartree-Fock을 구현할 때 이 행렬을 직접 만들어볼 것입니다.


## 다음 이야기

널리 쓰이는 근사법 중 하나인 변분법에 대해서 알아보았습니다. 다음 절에서는 또 다른 방법인 섭동론에 대해 살펴보겠습니다. 섭동론은 이미 정답을 알고 있는 문제에서 출발해서 다른 문제를 푸는 방법론입니다. 변분법과는 조금 다른 성격이죠. 하지만 마찬가지로 양자역학뿐만 아니라 다른 분야에서도 많이 쓰이는 방법입니다.


## 확인 문제
1. 위 코드에서 기저를 5개, 6개로 늘려보세요. 에너지가 계속 내려가나요? $\mathbf{S}$의 고유값도 함께 확인해보세요.
2. `eigh(H, S)`가 반환하는 두 번째와 세 번째 고유값은 무엇을 뜻할까요? 4장의 $E_2$, $E_3$와 비교해보세요.
3. 기저 함수를 규격화하지 않으면 어떻게 될까요? 코드에서 규격화 부분을 빼고 계산해보세요.
