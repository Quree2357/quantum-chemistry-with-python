# 15.1. Hückel 근사

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/15-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

1931년 독일의 물리학자 Hückel은 공액 구조를 가지는 분자들의 분자 오비탈을 계산할 수 있는 간단한 근사법을 제안했습니다. 이번 절에서는 이 방법이 어떻게 세워졌는지 알아보고 몇 가지 간단한 분자들에 적용해보겠습니다.


## Hückel 이론

Hückel 이론의 핵심 아이디어는 분자들의 화학적/분광학적 성질을 결정하는 $\pi$ 전자들만 따로 떼어서 계산한다는 것입니다. $\sigma$ 전자들은 분자의 골격을 구성하여 어떤 $\pi$ 오비탈들이 상호작용을 하는지만 결정하기 때문에 무시해도 된다는 것이죠. 이것을 $\sigma$-$\pi$ 분리라고 부릅니다.  

그리고 여기에 세 가지를 추가로 가정했습니다. 1) 모든 탄소의 Coulomb 적분은 같다. 2) 인접한 원자 사이의 공명 적분은 같은 값을 가지고, 인접하지 않으면 0이다. 3) 중첩 적분은 무시한다.  

이 방법은 분자를 이루는 원자들이 한 평면 위에 놓여 있을 때 성립하는 근사입니다. 분자 평면을 $xy$-평면이라고 두고, 분자 평면에 대한 거울 반사를 생각해보면 다른 오비탈들은 부호가 그대로겠지만 $p_z$ 오비탈은 부호가 바뀔 겁니다. 결국 $p_z$ 오비탈은 자기들끼리만 섞이고 $s$ 오비탈이나 $p_x$, $p_y$ 오비탈과는 섞이지 않겠죠. 그래서 Hückel 이론에서는 $\pi$ 오비탈만 따로 뗄 수 있습니다.

이렇게 가정하고 나면 Hamiltonian 행렬을 이렇게 쓸 수 있습니다.

$$
\mathbf{H} = \alpha \mathbf{I} + \beta \mathbf{A}
$$

여기서 $\mathbf{A}$는 인접 행렬(adjacency matrix)이라고 부르는데, 분자들의 결합이 어떻게 연결되어 있는지를 알려줍니다. 원자 $i$와 원자 $j$가 결합하고 있으면 $A_{rs}=1$이고 아니면 $0$인 모양이죠. 결국 원자들이 어떻게 연결되어 있는지만 알면 분자의 Hamiltonian 행렬이 정해지는 겁니다.


## 에틸렌과 부타디엔

가장 간단한 분자인 에틸렌부터 살펴봅시다. 공액 구조가 있는 건 아니지만 $\pi$ 결합을 가지는 가장 단순한 분자라서 시작점으로 삼기에는 충분합니다.

에틸렌 분자는 탄소 두 개가 이중 결합으로 연결되어 있는 구조입니다. 그러면 Hamiltonian 행렬은 다음과 같이 쓸 수 있죠.

$$
\mathbf{H} = \begin{bmatrix}
\alpha & \beta \\
\beta & \alpha
\end{bmatrix}
$$

영년 방정식 $(\alpha-E)^2 - \beta^2 = 0$을 풀면 에너지는 $\alpha \pm \beta$가 나옵니다. 에틸렌의 $\pi$ 전자는 두 개이니 모두 결합성 오비탈에 들어가겠군요. 따라서 $\pi$ 전자의 전체 에너지는 $2\alpha + 2\beta$가 됩니다. ($\beta$의 부호가 음수임에 주의하세요!)

탄소 네 개가 연결된 부타디엔은 어떨까요? 이 분자는 이중 결합과 단일 결합이 번갈아 가면서 존재하는 공액 구조를 가집니다.

$$
\mathbf{H} = \begin{bmatrix}
\alpha & \beta & 0 & 0 \\
\beta & \alpha & \beta & 0 \\
0 & \beta & \alpha & \beta \\
0 & 0 & \beta & \alpha
\end{bmatrix}
$$

마찬가지로 영년 방정식을 풀면 에너지가 나올 겁니다. 여기서 재미있는 건 방정식이 짝수 차수 항만으로 이루어져 있어서 부호가 반대인 쌍이 두 쌍 나온다는 것입니다.

$$
E_1 = \alpha + \frac{\sqrt{5} + 1}{2}\beta, \quad E_2 = \alpha + \frac{\sqrt{5} - 1}{2}\beta
$$
$$
E_3 = \alpha - \frac{\sqrt{5} + 1}{2}\beta, \quad E_4 = \alpha - \frac{\sqrt{5} - 1}{2}\beta
$$

부타디엔은 $\pi$ 전자가 4개이니 아래 두 오비탈에 2개씩 들어갈 겁니다. 따라서 $\pi$ 전자의 전체 에너지는 $4\alpha + 2\sqrt{5}\beta \approx 4\alpha + 4.472\beta$가 됩니다.

같은 방식으로 더 길게 연결된 분자들에도 적용할 수 있는데요. 코드로 만들어서 쉽게 계산할 수 있도록 해봅시다.

```python
import numpy as np


def huckel(n_carbons, bonds, alpha=0.0, beta=-1.0):
    H = alpha * np.eye(n_carbons)
    for i, j in bonds:
        H[i, j] = beta
        H[j, i] = beta
    E, C = np.linalg.eigh(H)
    return E, C


E, C = huckel(2, [(0, 1)])
print("에틸렌의 에너지 계수:", -E)

E, C = huckel(4, [(0, 1), (1, 2), (2, 3)])
print("부타디엔의 에너지 계수:", np.round(-E, 3))
```
```
에틸렌의 에너지 계수: [ 1. -1.]
부타디엔의 에너지 계수: [ 1.618  0.618 -0.618 -1.618]
```

만약 부타디엔 분자가 우리가 일반적으로 그리는 구조처럼 이중 결합과 단일 결합으로 이루어진 분자였다면 에너지가 $4\alpha + 4\beta$였을 겁니다. 하지만 실제로는 $0.472\beta$만큼 더 안정하죠. $\pi$ 전자가 결합 전체에 걸쳐서 퍼지면서 비편재화되기 때문에 생기는 현상입니다. 이 차이를 **비편재화 에너지(delocalization energy)**라고 부릅니다.  

각 분자 오비탈의 계수는 어떤지 볼까요?
```python
import matplotlib.pyplot as plt

C = C * np.sign(C[0])

fig, axes = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
for k, ax in enumerate(axes):
    colors = ["crimson" if c > 0 else "steelblue" for c in C[:, k]]
    ax.bar(range(1, 5), C[:, k], color=colors)
    ax.axhline(0, color="black", lw=1)
    ax.set_xticks(range(1, 5))
    ax.set_title(rf"$\alpha {-E[k]:+.3f}\beta$", fontsize=12)
    ax.set_xlabel("Carbon")
axes[0].set_ylabel("coefficient")
plt.tight_layout()
plt.show()
```
![부타디엔의 분자 오비탈 계수](/assets/image-106.png)

계수의 부호가 양수이면 $p_z$ 오비탈의 방향이 $+z$ 방향이라는 뜻이고, 음수이면 $-z$ 방향이라는 뜻입니다. 그림을 보면 알겠지만 마디의 개수가 0, 1, 2, 3개로 늘어나고 있고 계수의 크기도 무언가 사인 함수 비슷한 모양을 하고 있습니다. 어디서 많이 본 것 같지 않나요?  

네, 바로 상자 속 입자 모델에서 봤던 결과입니다! 이건 우연이 아니죠. 4장에서 시아닌 염료를 상자 속 입자 모델로 근사했던 것을 기억하시나요? 더 자세한 이야기는 다음 절에서 계속 해보도록 하죠.


## 전자가 얼마나 퍼져있을까

한 가지만 더 살펴보겠습니다. 위에서 봤던 분자 오비탈의 계수로부터 유용한 값 두 개를 계산할 수 있는데요. 바로 각 탄소 원자에 전자가 얼마나 존재하는지에 대한 전자 밀도와 결합이 얼마나 강한지를 나타내는 결합 차수입니다.

분자 오비탈 $k$에서 원자 $i$에 대한 계수를 $c_{ik}$라고 하고 그 오비탈에 있는 전자 수를 $n_k$라고 하면, $\pi$ 전자 밀도는 다음과 같이 계산할 수 있습니다.

$$
q_i = \sum_k n_k c_{ik}^2
$$

즉, 원자 $r$ 위에서 $\pi$ 전자를 평균적으로 몇 개 발견할 수 있는지를 나타내는 값이죠.

$\pi$ 결합 차수는 다음과 같이 계산합니다.

$$
p_{ij} = \sum_k n_k c_{ik} c_{jk}
$$

두 원자의 계수가 같은 부호이면 결합성이고 반대 부호이면 반결합성이 되겠죠. 그리고 당연히 전자가 들어있는 오비탈에 대해서만 합을 계산하면 됩니다.

부타디엔 분자에 대해서 두 값을 계산해봅시다.

```python
def pi_bonds(C, n_electrons):
    occ = np.zeros(C.shape[1])
    occ[: (n_electrons // 2)] = 2
    occ[n_electrons // 2] += n_electrons % 2
    q = (C**2) @ occ
    p = (C * occ) @ C.T
    return q, p


q, p = pi_bonds(C, 4)
print("pi 전자 밀도:", np.round(q, 3))
print("pi 결합 차수:", np.round(p[0, 1], 3), np.round(p[1, 2], 3), np.round(p[2, 3], 3))
```
```
pi 전자 밀도: [1. 1. 1. 1.]
pi 결합 차수: 0.894 0.447 0.894
```

$\pi$ 전자 밀도는 모든 탄소에서 1로 똑같습니다. 전자가 완전히 비편재화되어 있다는 뜻이죠. 탄소로만 이루어진 공액 구조에서는 항상 성립하는 성질입니다. 반면 $\pi$ 결합 차수는 조금 다르네요. 양쪽 말단의 두 결합은 0.894로 이중 결합에 가깝지만 가운데 결합은 0.447로 약 1.5중 결합에 가깝습니다. 이것도 비편재화 때문이죠.


## 다음 이야기

Hückel 이론을 에틸렌과 부타디엔과 같은 선형 분자에 적용해봤습니다. 길이가 좀 더 긴 분자들은 어떨까요? 그리고 벤젠 같은 고리형 공액 분자에서는 이 이론이 어떻게 적용될까요? 다음 절에서 4장과 7장에서 미뤄뒀던 이야기들을 정리해보겠습니다.


## 확인 문제
1. 에틸렌의 이중 결합의 길이는 1.34 Å, 에테인의 단일 결합의 길이는 1.54 Å입니다. 부타디엔의 세 결합의 결합 길이를 예측해보세요.
2. 알릴 라디칼(CH₂=CH-CH₂)은 탄소 세 개가 연결되어 있고 $\pi$ 전자가 3개입니다. `huckel` 함수로 에너지 준위를 구하고, 홀전자가 어느 탄소에 주로 있는지 HOMO의 계수로 확인해보세요.
3. 부타디엔의 비편재화 에너지는 $0.472 \beta$였습니다. 헥사트라이엔(탄소 6개짜리 사슬)에서는 얼마일까요? 사슬이 길어지면 어떻게 변할까요?
