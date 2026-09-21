# 13.2. 군의 표현과 기약표현

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/13-02.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

앞 절에서 분자의 대칭 조작과 대칭 요소에 대해 살펴봤습니다. 그리고 대칭 조작이 점군이라는 대수학적 구조를 이룬다는 것을 봤죠. 점군의 군 곱셈표를 만들 때는 행렬의 곱셈으로 계산했었습니다. 이렇게 군의 요소들을 행렬에 대응시키는 것을 **표현(representation)**이라고 합니다. 이번 절에서는 군의 표현에 대해서 조금 더 알아보도록 하죠.


## 군을 행렬로 '표현'한다

**군 표현(group representation)**이란 엄밀하게는 군과 벡터 공간의 선형 변환(linear transformation)을 대응시키는 것입니다. 수학적인 부분을 더 깊게 들어가지는 않겠습니다. 0.1.3에서 선형 변환과 행렬이 같은 것이라고 했었죠. 그러니 우리에게 지금 필요한 것은 행렬뿐입니다.  

대칭 조작은 일반적으로 교환 법칙이 성립하지 않는 것을 앞 절에서 암모니아 분자의 예시로 봤습니다. 행렬의 곱셈도 교환 법칙이 성립하지 않아서 점군을 행렬로 표현하는 것이 자연스럽죠. 이 경우를 특별하게 **행렬 표현(matrix representation)**이라고 합니다.  

여기서 주의할 점은, 대칭 조작을 행렬로 표현할 때 가능한 경우의 수가 무한히 많다는 것입니다. 똑같은 조작이라도 어느 방향에서 보냐에 따라 표현이 달라지거든요. 따라서 먼저 기저(basis)를 정해주는 것이 매우 중요합니다. 양자화학에서는 원자의 위치뿐만 아니라 오비탈에 대해서도 대칭 조작을 적용하기 때문에 기저를 벡터 대신 함수로 잡는 경우도 많습니다. 예를 들면 오비탈이나 진동 모드들의 집합과 같은 형태로요.  


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
s_N & s_B & s_C & s_A
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

*어, 보통 행렬이 열벡터 앞에 곱해지지 않아요?* 그렇긴 한데, 이 경우에는 관례적으로 기저를 행벡터 형태로 쓰고 변환 행렬을 뒤에 곱합니다. 이렇게 해야 행렬의 곱셈표와 대칭 조작의 곱셈표가 일치하거든요. 다시 말해, $D(R_1)D(R_2)=D(R_1R_2)$를 만족한다는 뜻입니다. 직접 확인해볼까요?
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
i = -np.eye(3)


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
13.1절의 확인 문제 4번의 답이기도 합니다. 잘 만들어졌나요? 이제 $s$ 오비탈 기저로 만든 표현 행렬도 확인해보죠.

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


## 행렬을 쪼개다

