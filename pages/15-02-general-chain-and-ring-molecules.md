# 15.2. 일반적인 사슬과 고리 분자

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/15-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>


## 폴리엔의 Hückel 근사

탄소 N개가 일렬로 연결되어 공액 구조를 이루는 폴리엔 분자를 생각해봅시다. 앞 절에서 봤던 에틸렌과 부타디엔의 Hückel 행렬의 모습을 생각해보면 일반화시킬 수 있겠죠.

$$
\mathbf{H} = \begin{bmatrix}
\alpha & \beta & 0 & & & & 0 \\
\beta & \alpha & \beta & & \cdots & & 0 \\
0 & \beta & \alpha & & & & 0 \\
 &  &  & \ddots & & & \vdots\\
 & \vdots & & & \alpha & \beta & 0 \\
 & & & & \beta & \alpha & \beta\\
0 & 0 & 0 & \cdots & 0 & \beta & \alpha
\end{bmatrix}
$$

0.2.3의 2차 차분 행렬처럼 삼중대각 행렬 형태로 나오게 됩니다. 그런데 행렬이 N차원이니 영년 방정식도 N차 방정식일텐데 이걸 대체 어떻게 풀어야 하죠?!  

다행히 삼중대각 행렬 중에서 특수한 모양을 가지는 행렬에 대한 영년 방정식의 해는 수학자들이 이미 찾아두었습니다. 더구나 우리의 Hamiltonian은 대칭 행렬이죠. 이 경우의 에너지와 계수는 다음과 같습니다.

$$
E_k = \alpha + 2\beta \cos \frac{k\pi}{N+1}
$$
$$
c_{ik} = \sqrt{\frac{2}{N+1}} \sin \frac{ik\pi}{N+1}
$$

흥미롭게도 삼각함수가 튀어나옵니다. 다시 한번 상자 속 입자가 생각나는 지점인데요. 직접 코드로 계산해서 상자 속 입자의 파동 함수와 겹쳐서 그려봅시다.

```python
import numpy as np
import matplotlib.pyplot as plt


def huckel(n_carbons, bonds, alpha=0.0, beta=-1.0):
    H = alpha * np.eye(n_carbons)
    for i, j in bonds:
        H[i, j] = beta
        H[j, i] = beta
    E, C = np.linalg.eigh(H)
    return E, C


def chain(N):
    return [(i, i + 1) for i in range(N - 1)]


N = 8
E, C = huckel(N, chain(N))
C = C * np.sign(C[0])

k = np.arange(1, N + 1)
i = np.arange(1, N + 1)

fig, axes = plt.subplots(1, 4, figsize=(12, 3), sharey=True)
x = np.linspace(0, N + 1, 100)
for k, ax in enumerate(axes):
    ax.bar(i, C[:, k], width=0.4, color="steelblue")
    ax.plot(x, np.sqrt(2 / (N + 1)) * np.sin((k + 1) * np.pi * x / (N + 1)), color="crimson", lw=2)
    ax.axhline(0, color="black", lw=1)
    ax.set_xlim(0, N + 1)
    ax.set_xticks(range(0, N + 2))
    ax.set_title(f"k = {k + 1}", fontsize=12)
    ax.set_xlabel("Carbon")
axes[0].set_ylabel("coefficient")
plt.tight_layout()
plt.show()
```
![폴리엔의 Hückel 근사](/assets/image-107.png)

파란색 막대로 나타난 것이 Hückel 계수, 곡선이 상자 속 입자의 파동 함수입니다. 이렇게 겹쳐서 보니 더욱 쉽게 알 수 있군요. Hückel 분자 오비탈은 바로 상자 속 입자의 파동 함수를 탄소 원자가 있는 자리에서만 찍어낸 것입니다.

그림에서 하나 더 볼 것이 있습니다. 상자 속 입자의 파동 함수가 0이 되는 상자의 양 끝점은 가장자리 탄소가 아니라 그 바깥쪽에 있습니다. 상자의 길이가 실제 결합의 총 길이보다 더 큰 거죠. 4.2절에서 상자의 길이를 정할 때 관례적으로 결합 하나씩을 더 준다고 했던 것 기억하시나요? Hückel 이론에서는 이 부분이 자연스럽게 나오게 됩니다.


## 시아닌과 폴리엔

4.2절에서 봤던 시아닌 염료에 Hückel 이론을 적용해보겠습니다. 실제로는 탄소 사슬의 양 끝에 질소 원자가 연결되어 있지만 일단은 모두 동일하다고 가정합니다. 아! 우리는 지금 $\beta$의 값을 모르네요... 적당히 3 eV로 정하고 계산해보죠.

```python
name = {1: "Cy3", 2: "Cy5", 3: "Cy7"}
obs = {1: 544, 2: 638, 3: 742}
box = {1: 453, 2: 579, 3: 706}

gaps = []
for k in [1, 2, 3]:
    n_atoms = 2 * k + 5
    n_pi = n_atoms + 1
    E, C = huckel(n_atoms, chain(n_atoms))
    homo = n_pi // 2 - 1
    gaps.append(E[homo + 1] - E[homo])

gaps = np.array(gaps)

print(f"{'':>3} {'gap':>6} {'Huckel':>7} {'box':>7} {'exp':>6}")
beta = 3
for k, g in zip([1, 2, 3], gaps):
    lam = 1239.84 / (g * beta)
    print(f"{name[k]:>3} {g:6.3f} {lam:7.1f} {box[k]:7.1f} {obs[k]:6d}")
```
```
       gap  Huckel     box    exp
Cy3  0.765   540.0   453.0    544
Cy5  0.618   668.7   579.0    638
Cy7  0.518   798.4   706.0    742
```

사실 실험값으로부터 $\beta$ 값을 맞춰서 썼습니다. 상자 속 입자 모델보다 훨씬 숫자가 잘 맞죠. 하지만 여기서 중요한 점은 두 모델이 같은 경향을 보인다는 것이고, 균일한 결합을 가지는 사슬에 Hückel 이론을 적용하면 상자 속 입자와 본질적으로 동일한 모델이라는 것입니다.  

4.2절의 TIP 박스에서 폴리엔은 결합 길이가 들쭉날쭉해서 상자 속 입자 모델을 쓰면 잘 안 맞는다고 했었죠. 그래서 폴리엔에 Hückel 이론을 적용할 때에는 $\beta$ 값을 두 개를 써야 합니다. 결합이 짧으면 오비탈이 더 많이 겹칠 테니 $\beta$ 값이 커지죠. 대충 1.1과 0.9 비율로 두고 사슬의 길이에 따른 HOMO-LUMO 간격을 비교해봅시다.

```python
def polyene(n, beta_double=-1.1, beta_single=-0.9):
    H = np.zeros((n, n))
    for i in range(n - 1):
        beta = beta_double if i % 2 == 0 else beta_single
        H[i, i + 1] = beta
        H[i + 1, i] = beta
    return np.linalg.eigvalsh(H)


ns = np.arange(2, 62, 2)
gap_equal, gap_alt = [], []
for n in ns:
    homo = n // 2 - 1
    E, _ = huckel(n, chain(n))
    gap_equal.append(E[homo + 1] - E[homo])
    E = polyene(n)
    gap_alt.append(E[homo + 1] - E[homo])

fig, ax = plt.subplots(figsize=(7, 5))
ax.plot(ns, gap_equal, "o--", color="steelblue", label="equal beta")
ax.plot(ns, gap_alt, "s--", color="crimson", label="alternating beta")
ax.axhline(0.4, color="crimson", ls=":", lw=2)
ax.set_xlabel("number of carbons")
ax.set_ylabel("HOMO-LUMO gap / |β|")
ax.set_ylim(0, 2.2)
ax.legend()
ax.grid(alpha=0.3)
plt.show()
```

$\beta$의 값을 하나만 썼을 때는 사슬의 길이가 길어질수록 HOMO와 LUMO 사이의 간격이 점점 줄어들어 0에 수렴합니다. 즉, 폴리엔의 길이가 길어질수록 흡수할 수 있는 빛의 에너지의 최소값이 점점 0에 가까워진다는 뜻입니다. 파장으로 보면 무한대로 증가하는 거죠.  

하지만 실제로는 탄소 사슬이 아무리 길어도 흡수 파장이 어느 이상으로 길어지지 않습니다. 빨간색 선으로 나타난 게 바로 이 경우죠. $\beta$의 값이 결합 길이에 따라서 번갈아 바뀌니 간격이 일정 값 이하로 떨어지지 않는 겁니다. 반면에 시아닌은 공명 구조 덕분에 결합 길이가 거의 일정해서 $\beta$ 값을 하나만 두어도 됐던 거고요.
