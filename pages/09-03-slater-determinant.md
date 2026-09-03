# 9.3. Slater 행렬식

앞 절에서 전자의 파동 함수는 반대칭이어야 한다고 했습니다. 그리고 전자가 두 개인 헬륨 원자에 대해서 반대칭인 파동 함수를 만들어봤죠. 그런데 리튬 원자는 전자가 세 개이니 항이 6개가 되고, 전자가 많아질수록 항의 개수가 급격하게 늘어납니다. 이걸 모두 손으로 일일이 쓸 수는 없습니다... 다전자 원자의 반대칭 파동 함수를 손쉽게 만들 수 있는 방법은 없을까요? 사실 우리가 이미 잘 알고 있는 개념을 이용하면 쉽게 만들 수 있습니다!


## 행렬식의 특징

선형대수학 시간에 행렬식(determinant)에 대해 배웠습니다. 그때는 역행렬을 계산하기 위한 수단 또는 연립방정식의 해에 대한 조건을 살펴보는 용도로 썼었죠. 행렬식을 어떻게 계산했는지 기억하시나요? 먼저 2 × 2 행렬의 경우를 살펴봅시다.

$$
\begin{vmatrix}
a & b \\
c & d
\end{vmatrix}
= ad - bc
$$

어라, 이렇게 보니까 앞 절에서 만들었던 반대칭 파동 함수와 모양이 비슷하네요! 실제로 이렇게 두면 정확히 같아집니다.

$$
\Psi = \frac{1}{\sqrt{2}}
\begin{vmatrix}
\psi_a(\mathbf{r}_1) & \psi_b(\mathbf{r}_1) \\
\psi_a(\mathbf{r}_2) & \psi_b(\mathbf{r}_2)
\end{vmatrix}
= \frac{1}{\sqrt{2}} \Big[ \psi_a(\mathbf{r}_1) \, \psi_b(\mathbf{r}_2) - \psi_a(\mathbf{r}_2) \, \psi_b(\mathbf{r}_1) \Big]
$$

이렇게 쓰는 방법을 **Slater 행렬식**이라고 부릅니다. 미국의 물리학자 Slater가 1929년에 도입한 방법이죠.  

행렬 안에 들어 있는 항을 잘 보면 다음과 같은 규칙이 있습니다.
- 같은 행에는 같은 전자 (1번 행에는 1번 전자)
- 같은 열에는 같은 스핀 오비탈 (1번 열에는 1번 스핀 오비탈)

선형대수학 시간에 배웠던 행렬식의 특징을 다시 떠올려보면 여러 가지가 있었는데 지금 우리에게 필요한 것은 두 개뿐입니다.  
첫 번째는 행렬의 두 행을 바꾸면 행렬식의 부호가 바뀐다는 것입니다. 행이 전자를 뜻하니 두 행을 바꾸는 것은 두 전자의 자리를 바꾸는 것이고 부호가 바뀐다는 것은 곧 반대칭성을 의미하는 것입니다. 그러니까 파동 함수를 Slater 행렬식으로 쓰면 파동 함수의 반대칭성이 저절로 보장되는 것이죠. 두 번째 성질은 행렬의 두 열이 같은 값을 가지면 행렬식이 0이 된다는 것입니다. 열은 스핀 오비탈을 뜻하니 두 열이 같다는 것은 두 전자가 같은 스핀 오비탈에 들어있다는 뜻입니다. 그러면 파동 함수가 0이 되죠. 행렬식의 성질로부터 우리가 원하는 파동 함수의 성질이 바로 얻어지는 것입니다.

직접 확인해볼까요?
```python
import sympy as sp

r1, r2 = sp.symbols("r1 r2")
psi_a, psi_b = sp.Function("psi_a"), sp.Function("psi_b")

M = sp.Matrix([[psi_a(r1), psi_b(r1)], [psi_a(r2), psi_b(r2)]])

print("행렬식:", sp.expand(M.det()))

M_swap = sp.Matrix([[psi_a(r2), psi_b(r2)], [psi_a(r1), psi_b(r1)]])

print("교환 후:", sp.expand(M_swap.det()))

M_same = sp.Matrix([[psi_a(r1), psi_a(r1)], [psi_a(r2), psi_a(r2)]])

print("같은 스핀 오비탈:", M_same.det())
```
```
행렬식: psi_a(r1)*psi_b(r2) - psi_a(r2)*psi_b(r1)
교환 후: -psi_a(r1)*psi_b(r2) + psi_a(r2)*psi_b(r1)
같은 스핀 오비탈: 0
```

Slater 행렬식이 파동 함수의 성질을 잘 표현하고 있습니다.


## 일반적인 Slater 행렬식

전자가 $N$개라면 어떻게 될까요? 여기서는 스핀 오비탈을 $\chi$라고 쓰겠습니다. (엄밀히 말하면 스핀 오비탈은 스핀 좌표에도 의존하므로 $\chi_1(\mathbf{x}_1)$처럼 쓰는 것이 정확하지만 여기서는 간단히 $\mathbf{r}$로 표기하겠습니다.)

$$
\Psi = \frac{1}{\sqrt{N!}}
\begin{vmatrix}
\chi_1(\mathbf{r}_1) & \chi_2(\mathbf{r}_1) & \cdots & \chi_N(\mathbf{r}_1) \\
\chi_1(\mathbf{r}_2) & \chi_2(\mathbf{r}_2) & \cdots & \chi_N(\mathbf{r}_2) \\
\vdots & \vdots & \ddots & \vdots \\
\chi_1(\mathbf{r}_N) & \chi_2(\mathbf{r}_N) & \cdots & \chi_N(\mathbf{r}_N)
\end{vmatrix}
$$

행렬 하나만으로 $N!$개의 항을 한 번에 표현할 수 있습니다. 이것이 Slater 행렬식 표현의 힘입니다.  
전자 3개인 경우를 직접 확인해봅시다.
```python
r3 = sp.symbols("r3")
chi_a, chi_b, chi_c = sp.Function("chi_a"), sp.Function("chi_b"), sp.Function("chi_c")

M3 = sp.Matrix(
    [[chi_a(r1), chi_b(r1), chi_c(r1)], [chi_a(r2), chi_b(r2), chi_c(r2)], [chi_a(r3), chi_b(r3), chi_c(r3)]]
)

slater = sp.expand(M3.det())

for term in slater.args:
    print(term)

print(f"\n항 개수: {len(slater.args)}")
```
```
chi_a(r1)*chi_b(r2)*chi_c(r3)
chi_a(r2)*chi_b(r3)*chi_c(r1)
chi_a(r3)*chi_b(r1)*chi_c(r2)
-chi_a(r1)*chi_b(r3)*chi_c(r2)
-chi_a(r2)*chi_b(r1)*chi_c(r3)
-chi_a(r3)*chi_b(r2)*chi_c(r1)

항 개수: 6
```
6개의 항이 부호까지 제대로 나왔습니다. 손으로 일일이 쓰려면 얼마나 귀찮을지 상상도 하기 싫네요.

[[TIP]]
행렬식 전체를 쓰는 것도 귀찮아 더 간략하게 만든 표기법도 있습니다. 어차피 전자도 $N$개고 스핀 오비탈도 $N$개니 그냥 대각선 성분만 쓰고 나머지는 알아서 생각하면 됩니다.

$$
\Psi = \ket{\chi_1,\chi_2,\cdots,\chi_N} = \ket{1,2,\cdots,N}
$$

Dirac의 브라켓 표기법에서 켓 기호만 따와서 썼죠. $\chi$를 쓰는 것도 귀찮아서 그냥 숫자만 쓰기도 합니다. 이 표기법은 17장 이후에 많이 보시게 될 겁니다. 예를 들어, 헬륨 원자의 바닥 상태는 $\ket{1s\alpha, 1s\beta}$ 또는 $\ket{1s, \bar{1s}}$라고 쓰게 됩니다.
[[/TIP]]


## 다음 이야기

Slater 행렬식 표현을 도입하면서 파동 함수의 반대칭성과 배타 원리가 자연스럽게 적용되고 $N!$개의 항을 한 번에 쓸 수 있게 되었습니다. 이 표현은 나중에 Hartree-Fock 이론을 세울 때 출발점으로 삼게 될 것입니다. 다만 아직 해결하지 못한 점도 있는데요, 오비탈 근사는 변수분리를 가능하게 만들기 위해서 억지로 곱 형태로 만든 거라 이 근사법을 사용하는 한 아직 정확도를 충분히 개선하지는 못합니다. 실제로는 전자들이 같은 부호의 전하를 띠기 때문에 서로를 피해 다니고 결국 전자들이 느끼는 퍼텐셜은 평균값이 아니라 순간적인 변화를 다 고려해주어야 합니다. 이것까지 담으려면 Slater 행렬식 하나만으로는 부족하고요. 얼마나 부족한지는 다음 절에서 확인해보겠습니다.


## 확인 문제
1. Slater 행렬식에서 두 열을 바꾸면 어떻게 될까요? 이것은 물리적으로 무슨 의미인가요?
2. 위 코드로 4 × 4 Slater 행렬식을 전개해보세요. 항이 몇 개가 나오나요?
3. 스핀오비탈 두 개가 정확히 같지는 않지만 아주 비슷하다면 행렬식이 어떻게 될까요? 이것이 무엇을 뜻할지 생각해보세요.
