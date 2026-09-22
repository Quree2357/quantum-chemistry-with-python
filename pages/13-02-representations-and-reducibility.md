# 13.2. 군의 표현과 기약표현

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/13-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞 절에서 분자의 대칭 조작과 대칭 요소에 대해 살펴봤습니다. 그리고 대칭 조작이 점군이라는 대수학적 구조를 이룬다는 것을 봤죠. 점군의 군 곱셈표를 만들 때는 행렬의 곱셈으로 계산했었습니다. 이렇게 군의 요소들을 행렬에 대응시키는 것을 **표현(representation)**이라고 합니다. 이번 절에서는 군의 표현에 대해서 조금 더 알아보도록 하죠.


## 군을 행렬로 '표현'한다

**군 표현(group representation)**이란 엄밀하게는 군과 벡터 공간의 선형 변환(linear transformation)을 대응시키는 것입니다. 수학적인 부분을 더 깊게 들어가지는 않겠습니다. 0.1.3에서 선형 변환과 행렬이 같은 것이라고 했었죠. 그러니 우리에게 지금 필요한 것은 행렬뿐입니다.  

대칭 조작은 일반적으로 교환 법칙이 성립하지 않는 것을 앞 절에서 암모니아 분자의 예시로 봤습니다. 행렬의 곱셈도 교환 법칙이 성립하지 않아서 점군을 행렬로 표현하는 것이 자연스럽죠. 이 경우를 특별하게 **행렬 표현(matrix representation)**이라고 합니다.  

여기서 주의할 점은, 대칭 조작을 행렬로 표현할 때 가능한 경우의 수가 무한히 많다는 것입니다. 똑같은 조작이라도 어느 방향에서 보느냐에 따라 표현이 달라지거든요. 따라서 먼저 기저(basis)를 정해주는 것이 매우 중요합니다. 양자화학에서는 원자의 위치뿐만 아니라 오비탈에 대해서도 대칭 조작을 적용하기 때문에 기저를 벡터 대신 함수로 잡는 경우도 많습니다. 예를 들면 오비탈이나 진동 모드들의 집합과 같은 형태로요.  


## 암모니아 분자의 대칭 조작의 행렬 표현

대칭 조작을 가하면 원자들의 위치가 바뀝니다. 원자의 좌표를 표기하려면 성분이 3개가 필요한데 물 분자만 해도 원자가 3개이니 3×3=9차원 벡터로 표시해야 하죠. 그런데 우리는 사실 분자 오비탈을 만들 때 어떤 원자 오비탈이 섞이는지 알고 싶었죠. 그러니까 오비탈 자체를 기저로 잡는 것이 편합니다.  

암모니아 분자의 경우를 예시로 들어보죠. 각 원자들의 $s$ 오비탈을 기저로 삼아서 4차원 벡터로 만들겠습니다. 그리고 대칭 조작을 했을 때 각 오비탈들이 어디로 가는지만 알면 되죠. 그러니 표현 행렬을 만드는 규칙이 아주 단순해집니다. $i$번째 오비탈이 $j$번째 오비탈 자리로 가면 $(j, i)$ 성분에 1을 쓰고, 나머지는 전부 0이 됩니다.

$$
[D(R)]_{ji} = \begin{cases}
1 & \text{오비탈 } i \text{가 오비탈 } j \text{의 자리로 갈 때}\\
0 & \text{그 외}
\end{cases}
$$

이렇게 자리를 바꾸기만 해주는 행렬을 치환 행렬(permutation matrix)이라고 부릅니다. 분자의 대칭 조작은 분자의 모양을 바꾸지 않으니 치환 행렬로 표현하는 게 훨씬 간단하고 편하죠.  

예를 들어, 회전 조작 $C_3$를 봅시다. 질소는 회전축 위에 있으니 가만히 있고, 수소 원자 세 개가 $A \to B \to C \to A$ 순서로 돌아가죠. 그러니까 이렇게 쓸 수 있습니다.

$$
\begin{bmatrix}
s_N & s_B & s_C & s_A
\end{bmatrix}
= \begin{bmatrix}
s_N & s_A & s_B & s_C
\end{bmatrix} D(C_3)
$$
$$
D(C_3) = \begin{bmatrix}
1 & 0 & 0 & 0 \\
0 & 0 & 0 & 1 \\
0 & 1 & 0 & 0 \\
0 & 0 & 1 & 0
\end{bmatrix}
$$

이렇게 두면 행렬의 곱셈표와 대칭 조작의 곱셈표가 일치합니다. 다시 말해, $D(R_1)D(R_2)=D(R_1R_2)$를 만족한다는 뜻입니다. 직접 확인해볼까요?
```python
import numpy as np


def C(n, axis="z"):
    """각 축에 대해 360/n도만큼 회전"""
    th = 2 * np.pi / n
    c, s = np.cos(th), np.sin(th)
    if axis == "z":
        return np.array([[c, -s, 0], [s, c, 0], [0, 0, 1.0]])
    if axis == "y":
        return np.array([[c, 0, s], [0, 1.0, 0], [-s, 0, c]])
    return np.array([[1.0, 0, 0], [0, c, -s], [0, s, c]])


def sigma_v(d):
    """z축을 포함하고 xy-평면과 각도 d를 이루는 평면에 대한 반사"""
    th = np.radians(d)
    return np.array([[np.cos(2 * th), np.sin(2 * th), 0], [np.sin(2 * th), -np.cos(2 * th), 0], [0, 0, 1.0]])


E = np.eye(3)

ops = {
    "E": E,
    "C3": C(3),
    "C3^2": C(3) @ C(3),
    "σ_a": sigma_v(0),
    "σ_b": sigma_v(120),
    "σ_c": sigma_v(240),
}


def which(M, table):
    for name, op in table.items():
        if np.allclose(M, op, atol=1e-9):
            return name
    return "??"


w = 7
print("대칭 조작의 군 곱셈표:")
print(" " * w + "".join(f"{k:>{w}}" for k in ops))
for a, A in ops.items():
    print(f"{a:>{w}}" + "".join(f"{which(A @ B, ops):>{w}}" for B in ops.values()))
```
```
대칭 조작의 군 곱셈표:
             E     C3   C3^2    σ_a    σ_b    σ_c
      E      E     C3   C3^2    σ_a    σ_b    σ_c
     C3     C3   C3^2      E    σ_c    σ_a    σ_b
   C3^2   C3^2      E     C3    σ_b    σ_c    σ_a
    σ_a    σ_a    σ_b    σ_c      E     C3   C3^2
    σ_b    σ_b    σ_c    σ_a   C3^2      E     C3
    σ_c    σ_c    σ_a    σ_b     C3   C3^2      E
```
```python
r, z = 0.9375, -0.3810
ammonia = np.array(
    [
        [0, 0, 0],
        [r, 0, z],
        [r * np.cos(2 * np.pi / 3), r * np.sin(2 * np.pi / 3), z],
        [r * np.cos(4 * np.pi / 3), r * np.sin(4 * np.pi / 3), z],
    ]
)
labels = ["N", "H", "H", "H"]
basis = ["s_N", "s_A", "s_B", "s_C"]


def orbital_rep(op, coords, labels):
    """s 오비탈 기저에 대한 표현 행렬"""
    n = len(coords)
    M = np.zeros((n, n))
    moved = coords @ op.T
    for i, (p, l) in enumerate(zip(moved, labels)):
        for j, (q, m) in enumerate(zip(coords, labels)):
            if l == m and np.allclose(p, q, atol=1e-6):
                M[j, i] = 1
                break
    return M


orb_reps = {k: orbital_rep(v, ammonia, labels) for k, v in ops.items()}
print("표현 행렬의 군 곱셈표:")
print(" " * w + "".join(f"{k:>{w}}" for k in orb_reps))
for a, A in orb_reps.items():
    print(f"{a:>{w}}" + "".join(f"{which(A @ B, orb_reps):>{w}}" for B in orb_reps.values()))
```
```
표현 행렬의 군 곱셈표:
             E     C3   C3^2    σ_a    σ_b    σ_c
      E      E     C3   C3^2    σ_a    σ_b    σ_c
     C3     C3   C3^2      E    σ_c    σ_a    σ_b
   C3^2   C3^2      E     C3    σ_b    σ_c    σ_a
    σ_a    σ_a    σ_b    σ_c      E     C3   C3^2
    σ_b    σ_b    σ_c    σ_a   C3^2      E     C3
    σ_c    σ_c    σ_a    σ_b     C3   C3^2      E
```
두 경우의 곱셈표가 완벽하게 똑같습니다. 행렬의 크기가 다른데도 같은 군을 표현할 수 있는 것이죠.


## 터져라, 시냅스! 갈라져라, 행렬!

물 분자나 암모니아 분자처럼 크기가 작으면 상관이 없지만 분자가 커지면 원자의 갯수만큼 행렬의 차원이 커집니다. 들고 다니기 편하게 쪼갤 수 있는 방법은 없을까요?  

암모니아 분자를 다시 봅시다. 위에서처럼 놓으면 회전 조작과 반사 조작 모두에 대해서 질소 원자의 위치는 변하지 않죠. 회전축과 거울면이 모두 질소 원자를 지나가니까요. 그러니 $s_N$ 오비탈은 다른 오비탈과 섞이지 않습니다. 그러니까 표현 행렬에서 질소는 따로 1×1 행렬로 떼어낼 수 있습니다.  

수소는 어떤가요? 회전 조작을 하면 서로 자리를 바꾸니 이대로는 더 쪼갤 수 없습니다. 하지만 만약 회전 조작에 대해서 바뀌지 않는 기저를 하나라도 가지도록 만들면 어떻게 될까요? 다음 기저를 생각해봅시다.

$$
(s_A + s_B + s_C, \, 2s_A - s_B - s_C, \, s_B - s_C)
$$

세 함수는 모두 선형독립이니 기저를 이룰 수 있습니다. 여기서 눈에 띄는 것은 첫 번째 기저가 모든 대칭 조작에 대해 불변이라는 것입니다. 그러니 이 기저도 따로 떼어낼 수 있겠군요. 이제 남은 것은 2차원입니다. 이런 기저를 **대칭 적응 기저(symmetry-adapted basis)**라고 부릅니다.  

이 기저에서 대칭 조작의 표현 행렬이 어떻게 나타나는지 살펴봅시다. 기저가 다르면 행렬 표현이 달라지니까요.
```python
# 규격화된 새로운 기저
new_basis = [
    ("sN", np.array([1, 0, 0, 0])),
    ("sA + sB + sC", np.array([0, 1 / np.sqrt(3), 1 / np.sqrt(3), 1 / np.sqrt(3)])),
    ("2sA - sB - sC", np.array([0, 2 / np.sqrt(6), -1 / np.sqrt(6), -1 / np.sqrt(6)])),
    ("sB - sC", np.array([0, 0, 1 / np.sqrt(2), -1 / np.sqrt(2)])),
]

new_basis_names = [b[0] for b in new_basis]
T = np.column_stack([b[1] for b in new_basis])

print("T는 직교 행렬인가:", np.allclose(T.T @ T, np.eye(4)))
print()
w = 16
for r in ["C3", "σ_a"]:
    print(f"새 기저에서 {r}의 표현 행렬")
    print(" " * w + "".join(f"{k:>{w}}" for k in new_basis_names))
    for b, row in zip(new_basis_names, T.T @ orb_reps[r] @ T):
        clean = [0.0 if abs(v) < 1e-12 else v for v in row]
        print(f"{b:>{w}}" + "".join(f"{v:>16.3f}" for v in clean))
    print()
```
```
T는 직교 행렬인가: True

새 기저에서 C3의 표현 행렬
                              sN    sA + sB + sC   2sA - sB - sC         sB - sC
              sN           1.000           0.000           0.000           0.000
    sA + sB + sC           0.000           1.000           0.000           0.000
   2sA - sB - sC           0.000           0.000          -0.500          -0.866
         sB - sC           0.000           0.000           0.866          -0.500

새 기저에서 σ_a의 표현 행렬
                              sN    sA + sB + sC   2sA - sB - sC         sB - sC
              sN           1.000           0.000           0.000           0.000
    sA + sB + sC           0.000           1.000           0.000           0.000
   2sA - sB - sC           0.000           0.000           1.000           0.000
         sB - sC           0.000           0.000           0.000          -1.000
```

표현 행렬이 같은 모양으로 나타나고 있습니다. 왼쪽 위부터 1×1 행렬 두 개와 2×2 행렬 하나가 대각선 방향으로 늘어서 있죠. 이런 모양의 행렬을 블록 대각 행렬(block diagonal matrix)이라고 부릅니다. 표기로는 이렇게 쓰기도 합니다.

$$
D^{(4)}=\underbrace{D^{(1)}}_{s_N} \oplus
\underbrace{D^{(1)}}_{s_A+s_B+s_C} \oplus
\underbrace{D^{(2)}}_{나머지}
$$

마지막 행렬의 원소들은 대칭 조작에 따라 달라지기는 하지만 우리는 이제 성분 16개 대신 6개만 들고 다니면 됩니다.


## 행렬, 어디까지 쪼개지니?

암모니아 분자의 표현 행렬은 기저를 잘 잡아서 (1, 1, 2)차원 조각들로 쪼갤 수 있었습니다. 만약 다른 기저를 택한다면 더 작은 조각으로 쪼갤 수 있을까요? 답은 '아니다'입니다. 마지막에 나타난 2차원 행렬은 어떤 기저를 잡아도 더 이상 쪼갤 수 없거든요. 이렇게 더 쪼갤 수 없는 경우를 **기약 표현(irreducible representation)**이라고 합니다.  

기약 표현은 점군의 고유한 성질입니다. 같은 점군은 같은 기약 표현을 갖죠. 그래서 이 기약 표현을 이용해 점군들을 또 분류할 수 있고, 이렇게 분류된 것을 **대칭 종(symmetry species)**이라고 합니다. 다음 절에서 대칭 종의 분류에 대해 살펴보고 이것이 왜 중요한지 알아보겠습니다.


## 다음 이야기

암모니아 분자의 표현 행렬을 잘 쪼개서 성분 6개만 들고 다녀도 되도록 했죠. 그런데 사실 이것도 분자가 작아서 다행이지 조금만 커져도 섞이는 오비탈들의 수가 많아져서 들고 다녀야 하는 양이 많아집니다. 그러니 조금 더 편하게 쓸 수 있었으면 좋겠습니다. 조각 하나를 숫자 하나로 부를 수 있다면 얼마나 좋을까요? 사실 그게 가능하고, 다음 절에서 살펴볼 내용입니다.


## 확인 문제

1. $σ_b$와 $C_3^2$ 조작의 4×4 표현 행렬을 만들어보세요. 어느 오비탈이 제자리에 있고 어느 오비탈이 자리를 바꾸나요?
2. 대칭 조작의 표현 행렬들의 대각합을 보세요. 13.1절에서 봤던 같은 유형에 속한 조작들의 값은 어떤가요?
3. $2s_A - s_B - s_C$ 대신 $2s_B - s_A - s_C$를 써도 될까요? 직교성을 확인해보고 행렬이 쪼개지는 것을 코드로 확인해보세요.
4. 물 분자의 표현 행렬을 만들어보세요. 암모니아의 경우와 달리 1×1 행렬만으로 완전히 쪼개지는데, 그 이유를 $C_{2v}$가 가환군이라는 사실과 연관지어 설명해보세요.
