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


## 강강수월래 를! 하는 탄소
 
이번에는 탄소 $N$개를 고리로 이어봅시다. 사슬과 달리 고리에는 끝이 없으니 7.1절의 고리 위 입자처럼 주기적 경계 조건이 됩니다. 이 경우에도 고유값을 식으로 쓸 수 있습니다.
 
$$
E_m = \alpha + 2\beta \cos\frac{2\pi m}{N} \qquad (m = 0, \pm1, \pm2, \dots)
$$
 
$m$과 $-m$이 같은 에너지를 주니 $m=0$을 빼면 준위가 두 개씩 축퇴되어 있습니다. 7.1절의 고리 위 입자와 똑같죠. 그리고 $m$이 작을 때 코사인을 전개하면 $E_m \approx \alpha + 2\beta - \beta\left(\frac{2\pi}{N}\right)^2 m^2$이 되어, 고리 위 입자처럼 에너지가 $m^2$에 비례합니다.
 
벤젠($N=6$)에 대해 계산하고, 13.5절에서 사영 연산자로 만든 SALC를 기저로 Hamiltonian 행렬을 다시 써봅시다. 13.3절에서 본 것처럼 기저를 바꾸는 것은 $\mathbf{U}^T\mathbf{H}\mathbf{U}$를 계산하는 것이었죠.
```python
def ring(n):
    return chain(n) + [(n - 1, 0)]
 
 
E, C = huckel(6, ring(6))
print("벤젠 x:", np.round(-E, 3))
 
# 13.5절에서 사영 연산자로 만든 SALC (A1, E1, E1, E2, E2, B1)
salc = np.array([
    [1, 1, 1, 1, 1, 1],
    [2, 1, -1, -2, -1, 1],
    [0, 1, 1, 0, -1, -1],
    [2, -1, -1, 2, -1, -1],
    [0, 1, -1, 0, 1, -1],
    [1, -1, 1, -1, 1, -1],
], dtype=float)
U = (salc / np.linalg.norm(salc, axis=1, keepdims=True)).T
 
H = np.zeros((6, 6))
for i, j in ring(6):
    H[i, j] = H[j, i] = -1.0
print("SALC 기저의 Hamiltonian (x로 표시):")
print(np.round(-(U.T @ H @ U), 3) + 0.0)
```
```
벤젠 x: [ 2.  1.  1. -1. -1. -2.]
SALC 기저의 Hamiltonian (x로 표시):
[[ 2.  0.  0.  0.  0.  0.]
 [ 0.  1.  0.  0.  0.  0.]
 [ 0.  0.  1.  0.  0.  0.]
 [ 0.  0.  0. -1.  0.  0.]
 [ 0.  0.  0.  0. -1.  0.]
 [ 0.  0.  0.  0.  0. -2.]]
```
 
벤젠의 에너지는 $\alpha + 2\beta$, $\alpha + \beta$ (2개), $\alpha - \beta$ (2개), $\alpha - 2\beta$입니다. 그리고 SALC를 기저로 쓰면 Hamiltonian 행렬이 대각 성분만 남습니다. 대각화를 하지 않았는데도 이미 대각 행렬인 거죠.
 
13.4절에서 대칭 종이 다른 오비탈끼리는 섞이지 않는다고 했습니다. 그런데 벤젠의 SALC 여섯 개 중에서 같은 대칭 종을 가지는 것은 $E_1$ 두 개, $E_2$ 두 개뿐이고, 이 둘은 각각 축퇴되어 있으니 서로 섞여도 에너지가 같습니다. 결국 행렬에 남을 수 있는 것은 대각 성분뿐입니다. 14.4절의 물 분자에서는 $A_1$ 오비탈이 세 개라서 3×3 블록을 따로 풀어야 했던 것과 대조적이죠. 대칭성만으로 벤젠의 분자 오비탈을 완전히 결정한 셈입니다.
 
13.5절에서 벤젠의 SALC 계수가 $\cos m\phi_k$, $\sin m\phi_k$ 모양이었고, $A_1$, $E_1$, $E_2$, $B_1$이 각각 $m = 0, \pm1, \pm2, 3$에 해당한다고 했었습니다. 이제 에너지도 붙일 수 있습니다. 위 식에 $N=6$을 넣으면 $m=0$은 $2\beta$, $m=\pm1$은 $\beta$, $m=\pm2$는 $-\beta$, $m=3$은 $-2\beta$가 나오죠. 고리 위 입자의 파동 함수를 탄소 자리에서 찍은 SALC가 그대로 분자 오비탈이고, 에너지는 $m$이 커질수록 올라갑니다.
 
[[TIP]]
13.5절에서는 부호가 바뀌지 않는 조작만 모은 $C_{6v}$를 썼지만, 벤젠의 실제 점군은 $D_{6h}$입니다. $D_{6h}$의 이름으로 쓰면 벤젠의 $\pi$ 분자 오비탈은 에너지 순서대로 $a_{2u}$, $e_{1g}$, $e_{2u}$, $b_{2g}$입니다. 다른 교재나 논문에서는 이 이름으로 더 자주 만날 겁니다.
[[/TIP]]
 
 
## 방향족성과 $4n+2$
 
벤젠의 $\pi$ 전자 여섯 개는 $\alpha + 2\beta$에 두 개, $\alpha + \beta$에 네 개가 들어가서 $\pi$ 전자의 전체 에너지가 $6\alpha + 8\beta$입니다. 에틸렌 세 개($6\alpha + 6\beta$)보다 $2\beta$만큼 더 안정하죠. 15.1절에서 본 부타디엔의 비편재화 에너지가 $0.472\beta$였으니 훨씬 큽니다. 이 큰 안정화가 벤젠을 특별하게 만드는, 이른바 **방향족성(aromaticity)**의 정체입니다.
 
7.1절에서 고리 위 입자에 전자를 채우면 2, 6, 10, … 개에서 껍질이 닫힌다고 했었습니다. $m=0$에 두 개, 그다음부터는 축퇴된 준위 하나에 네 개씩 들어가니 $4n+2$죠. Hückel 고리도 축퇴 구조가 똑같으니 같은 규칙이 나옵니다. 그러면 전자 수가 $4n+2$가 아닐 때는 어떻게 될까요? 고리 크기와 전자 수를 바꿔가며 홀전자 수를 세어봅시다.
```python
def fill(x, n_pi):
    occ = np.zeros(len(x))
    filled = 0
    for level in sorted(set(np.round(x, 6)), reverse=True):
        idx = np.where(np.isclose(x, level))[0]
        k = min(2 * len(idx), n_pi - filled)
        occ[idx] = k / len(idx)
        filled += k
    return occ
 
 
molecules = [
    ("C4H4", 4, 4),
    ("C5H5-", 5, 6),
    ("C6H6", 6, 6),
    ("C7H7+", 7, 6),
    ("C8H8", 8, 8),
]
 
print(f"{'':>6} {'ring':>4} {'pi e':>5} {'unpaired':>9}")
for name_m, n, n_pi in molecules:
    E, _ = huckel(n, ring(n))
    occ = fill(-E, n_pi)
    unpaired = sum(min(o, 2 - o) for o in occ)
    print(f"{name_m:>6} {n:4d} {n_pi:5d} {unpaired:9.0f}")
```
```
       ring  pi e  unpaired
  C4H4    4     4         2
 C5H5-    5     6         0
  C6H6    6     6         0
 C7H7+    7     6         0
  C8H8    8     8         2
```
 
$\pi$ 전자가 6개인 $\text{C}_5\text{H}_5^-$, $\text{C}_6\text{H}_6$, $\text{C}_7\text{H}_7^+$는 고리의 크기와 상관없이 홀전자가 없습니다. 반면 $\pi$ 전자가 4개인 사이클로뷰타다이엔과 8개인 사이클로옥타테트라엔은 홀전자가 두 개씩 생기죠.
 
사이클로뷰타다이엔을 자세히 봅시다. 에너지 준위는 $\alpha + 2\beta$, $\alpha$ (2개), $\alpha - 2\beta$입니다. 전자 네 개 중 두 개는 $\alpha + 2\beta$에 들어가고, 나머지 두 개는 축퇴된 $\alpha$ 준위 두 개에 하나씩 들어갑니다. 14.2절의 산소 분자처럼 Hund 규칙에 따라서요. $\pi$ 전자의 전체 에너지는 $4\alpha + 4\beta$로 에틸렌 두 개와 똑같으니 비편재화로 얻는 안정화도 전혀 없습니다. 오히려 축퇴된 준위가 반만 채워진 불안정한 상태가 되죠. 이런 분자를 **반방향족(antiaromatic)**이라고 부릅니다. 실제로 사이클로뷰타다이엔은 극저온에서만 겨우 만들 수 있을 정도로 불안정하고, 사이클로옥타테트라엔은 평면 구조를 포기하고 욕조 모양으로 구부러져서 이 문제를 피합니다.
 
[[TIP]]
7.1절에서 벤젠을 고리 위 입자로 계산했을 때 흡수 파장은 212.6 nm였고 실제로는 255 nm였습니다. Hückel 이론에서 벤젠의 HOMO–LUMO 간격은 $2|\beta|$인데, 시아닌에서 맞춘 $|\beta| = 3.08$ eV를 넣으면 약 201 nm가 나옵니다. 오히려 더 멀어졌죠. 사실 벤젠의 255 nm 흡수는 HOMO($E_1$)에서 LUMO($E_2$)로 가는 하나의 전이가 전자 사이의 반발 때문에 에너지가 다른 여러 상태로 갈라진 것 중 가장 낮은 것입니다. 실제로 벤젠은 약 200 nm와 180 nm에서도 빛을 흡수하죠. 전자 반발을 아예 무시하는 Hückel 이론으로는 이 갈라짐을 설명할 수 없습니다. 이 이야기는 24장에서 다시 하겠습니다.
[[/TIP]]
 
 
## 다음 이야기
 
지금까지는 탄소로만 이루어진 분자를 다뤘습니다. 다음 절에서는 질소나 산소 같은 다른 원자가 들어가면 무엇이 달라지는지 살펴보겠습니다. 14.3절에서 미뤄뒀던 push-pull 염료 이야기도 여기서 할 수 있습니다. 그리고 Hückel 이론이 어디까지 믿을 만한지도 정리해보겠습니다.
 
 
## 확인 문제
 
1. 탄소 8개 사슬에 대해 위의 에너지 식으로 HOMO–LUMO 간격을 구하고, 상자 속 입자에서 $L = 9d$로 둔 결과와 비교해보세요. 사슬이 아주 길어지면 두 결과는 어떻게 될까요?
2. `polyene` 함수에서 `b_double`과 `b_single`을 바꿔가며 긴 사슬의 HOMO–LUMO 간격이 $2|\beta_1 - \beta_2|$로 수렴하는지 확인해보세요.
3. 나프탈렌은 벤젠 고리 두 개가 한 변을 공유한 분자입니다. 탄소 10개의 연결 정보를 직접 적어서 `huckel` 함수로 $\pi$ 에너지 준위를 구하고, 비편재화 에너지를 벤젠 두 개와 비교해보세요.