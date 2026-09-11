# 11.4. 헬륨 원자의 근사

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/11-04.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

9장에서 헬륨 원자의 Schrödinger 방정식을 보면서 정확한 해석해를 구하는 방법은 없다고 했었죠. 오비탈 근사를 쓰면 오차가 거의 30 eV 나서 쓸 수가 없는 값이 나왔었고요. 그리고 10장에서 가리움과 유효 핵 전하라는 개념을 도입하면서 전자들이 느끼는 평균적인 퍼텐셜 에너지로 따져보면 조금 더 실험값들과 잘 맞았었죠. 그런데 그게 대체 왜 그렇게 나오는지는 알지 못했습니다. 헬륨 원자의 유효 핵 전하가 왜 1.70인지 같은 것들이요. 이번 절에서는 방금 보았던 변분법과 섭동론의 두 가지 근사법을 통해 헬륨 원자의 바닥 상태 에너지를 계산해보고 유효 핵 전하의 정체도 밝혀보겠습니다.


## 헬륨 원자와 섭동

헬륨 원자의 Hamiltonian을 다시 봅시다. 원자 단위계로 썼었죠.

$$
\hat{H} = \left( -\frac{1}{2} \nabla_1^2 - \frac{2}{r_1} \right) + \left( -\frac{1}{2} \nabla_2^2 - \frac{2}{r_2} \right) + \frac{1}{r_{12}}
$$

앞의 항들은 수소꼴 원자 두 개에 대한 Hamiltonian이니 정확히 풀 수 있는데 마지막에 있는 전자 간 상호작용 항 때문에 변수분리를 못한다고 했었습니다. 그런데 섭동론의 관점에서 보면 지금 이 형태가 딱 맞습니다. 마지막 항을 섭동으로 놓으면 되죠.

$$
\hat{H}^{(0)} = \left( -\frac{1}{2} \nabla_1^2 - \frac{2}{r_1} \right) + \left( -\frac{1}{2} \nabla_2^2 - \frac{2}{r_2} \right) \qquad \hat{H}'= \frac{1}{r_{12}}
$$

비섭동 Hamiltonian에 대한 해는 이미 알고 있습니다. 9장에서 오비탈 근사를 했을 때와 똑같은 결과입니다.

$$
\psi^{(0)} = \sqrt{\frac{8}{\pi}}e^{-2r_1} \cdot \sqrt{\frac{8}{\pi}}e^{-2r_2} \qquad E^{(0)} = -4 \left( \frac{1}{2n_1^2} + \frac{1}{2n_2^2} \right)
$$

두 전자가 모두 $1s$ 오비탈에 있다고 하면 $n=1$이니까 $E^{(0)} = -4$ Hartree가 되죠. 9.1절에서 계산했던 -108.85 eV와 같은 값입니다. 이제 섭동에 대한 보정 항을 계산하면 되겠네요.

에너지의 1차 보정 항은 다음과 같습니다.

$$
E^{(1)} = \left\langle \psi^{(0)} \middle| \hat{H}' \middle| \psi^{(0)} \right\rangle = \left\langle \psi^{(0)} \middle| \frac{1}{r_{12}} \middle| \psi^{(0)} \right\rangle
$$

이 적분은 좀 까다롭습니다. 두 전자의 좌표 각각에 대해서 적분을 해야 하니 6중 적분이 나오고 거리 항이 분모에 있거든요. 계산 과정은 생략하고 결과만 보면 $E^{(1)} = \frac{5}{4}$ Hartree가 나옵니다. 보정하지 않은 에너지에 더해주면 -74.83 eV를 얻습니다. 실험값이 -79.00 eV였던 것을 생각하면 많이 가까워지기는 했지만 아직도 오차가 조금 있네요. 반발 항의 기여분이 34 eV나 되니 전체의 30%가 넘습니다. 섭동이 이렇게 크면 오차가 잘 줄어들지 않습니다.

그러면 이제 2차 보정 항까지 고려해봐야겠네요. 1차 보정을 더해주었을 때 실험값보다 오히려 높아져버렸으니 음수로 나오는 2차 보정을 더해주면 다시 낮아질 것입니다. 다만 이 계산은 정말 어렵습니다. 중간 과정을 쓰는 것조차 매우 복잡하니 문헌에서 계산되어 있는 값을 가져오겠습니다. $E^{(2)} = -0.15766$ Hartree입니다. 이 값을 더해주면 2차 보정 에너지는 -79.12 eV가 됩니다. 오차가 0.12 eV까지 줄어들었네요! 그런데 이번에는 실험값보다 작아졌습니다. 11.3절에서 얘기했던 것처럼 섭동론에는 상한이나 하한에 대한 보장이 없다는 사실을 봤습니다.

[[TIP]]
3차 보정 항까지 고려해주면 에너지가 -79.00 eV가 됩니다. 거의 실험값과 일치하죠. 13차 보정까지 계산한 문헌도 있는데 소수점 아래 6번째 자리까지 일치하는 결과를 보여줬습니다. 이 문헌에서는 변분법과 섭동론을 섞은 방법을 사용했는데 무려 100개의 항을 가진 시험 함수를 썼다고 하네요. (C. W. Scherr와 R. E. Knight의 논문 "Two-Electron Atoms III. A Sixth-Order Perturbation Study of the 1⁢^1s Ground State"을 참고하세요.) 이렇게 적절한 근사법을 사용하면 매우 정확한 수준의 값을 얻을 수 있다는 것을 보여줍니다.
[[/TIP]]


## 헬륨 원자와 변분법

이번에는 변분법 관점에서 접근해보겠습니다. 10.1절에서 봤던 것처럼 가리움 효과 때문에 전자가 실제로 느끼는 유효 핵 전하는 2보다 작습니다. 그러면 핵이 전자를 당기는 인력이 줄어드니 궤도가 조금 더 부풀어야 하겠죠. 그러니 비섭동 항에서 잡았던 것처럼 $Z=2$로 두는 대신 이 값을 매개변수로 두는 것이 자연스럽습니다.

$$
\psi_t = \underbrace{\sqrt{\frac{Z^3}{\pi}}e^{-Zr_1}}_{\psi_1} \cdot \underbrace{\sqrt{\frac{Z^3}{\pi}}e^{-Zr_2}}_{\psi_2} = \frac{Z^3}{\pi}e^{-Z(r_1 + r_2)}
$$

이 함수를 시험 함수로 두고 변분 원리를 적용하면 되죠. 적분을 좀 해야 하는데 이번에는 SymPy에게 시켜봅시다. 일단 분모부터 계산해보죠.

$$
\int \psi_t^* \psi_t \, d\tau = \int \psi_1^* \psi_1 \, d\tau \int \psi_2^* \psi_2 \, d\tau
$$

```python
import sympy as sp

Z, r, r1, r2 = sp.symbols("Z r r_1 r_2", positive=True)

phi_1 = sp.sqrt(Z**3 / sp.pi) * sp.exp(-Z * r1)
phi_2 = sp.sqrt(Z**3 / sp.pi) * sp.exp(-Z * r2)

integral1 = sp.integrate(phi_1.conjugate() * phi_1 * 4 * sp.pi * r1**2, (r1, 0, sp.oo))
integral2 = sp.integrate(phi_2.conjugate() * phi_2 * 4 * sp.pi * r2**2, (r2, 0, sp.oo))
print("규격화 확인:", integral1 * integral2)
```
```
규격화 확인: 1
```
오, 규격화가 되어 있습니다! 분모는 이제 무시해도 되겠네요. 분자를 살펴봅시다. 먼저 Hamiltonian이 어떻게 생겼는지 다시 보죠.

$$
\hat{H} = \left( -\frac{1}{2} \nabla_1^2 - \frac{2}{r_1} \right) + \left( -\frac{1}{2} \nabla_2^2 - \frac{2}{r_2} \right) + \frac{1}{r_{12}}
$$

이걸 가운데 넣고 적분해야 한다니 생각보다 끔찍한데요... 이렇게 써보면 어떨까요?

$$
\hat{H} = \left( -\frac{1}{2} \nabla_1^2 - \frac{Z}{r_1} \right) + \left( -\frac{1}{2} \nabla_2^2 - \frac{Z}{r_2} \right) + \frac{Z-2}{r_1} + \frac{Z-2}{r_2} + \frac{1}{r_{12}}
$$

항이 조금 더 많아지긴 했는데 앞의 두 항이 각각 핵 전하가 $Z$인 수소꼴 원자의 Hamiltonian이 됐죠. 그리고 우리가 정한 시험 함수는 딱 이 경우에 맞는 파동 함수 모양입니다. 그러니 Hamiltonian을 취하면 고유값이 나오겠죠.

$$
\left[ \left( -\frac{1}{2} \nabla_1^2 - \frac{Z}{r_1} \right) + \left( -\frac{1}{2} \nabla_2^2 - \frac{Z}{r_2} \right) \right] \psi_t = -Z^2 \psi_t
$$

그러니까 나머지는 이렇게 쓸 수 있습니다.

$$
\int \psi_t^* \hat{H} \psi_t \, d\tau = -Z^2 \int \psi_t^* \psi_t \, d\tau + (Z-2)\left[\int \frac{\psi_t^* \psi_t}{r_1} \, d\tau + \int \frac{\psi_t^* \psi_t}{r_2} \, d\tau \right]+ \int \frac{\psi_t^* \psi_t}{r_{12}} \, d\tau
$$

이제 각 항의 적분을 해봅시다. 맨 처음 적분은 분모와 같은 모양인데 아까 1이라고 했었죠. 그 다음 적분부터만 계산하면 됩니다.

```python
integral3 = sp.integrate(phi_1 * phi_1 / r1 * 4 * sp.pi * r1**2, (r1, 0, sp.oo))
integral4 = sp.integrate(phi_2 * phi_2 / r2 * 4 * sp.pi * r2**2, (r2, 0, sp.oo))
print("첫 번째 적분값: ", integral3 * integral2)
print("두 번째 적분값: ", integral4 * integral1)
```
```
첫 번째 적분값:  Z
두 번째 적분값:  Z
```
마지막 항이 조금 문제입니다. $r_{12}$가 $r_1$과 $r_2$에 모두 의존하니 SymPy로 한 번에 쓰기가 어렵거든요. 그래서 경우를 둘로 나눠서 계산합니다.

```python
f = (4 * sp.pi) ** 2 * phi_1**2 * phi_2**2 * r1**2 * r2**2
case_1 = sp.integrate(f / r1, (r2, 0, r1))
case_2 = sp.integrate(f / r2, (r2, r1, sp.oo))
integral5 = sp.integrate(sp.simplify(case_1 + case_2), (r1, 0, sp.oo))
print("세 번째 적분값: ", integral5)
```
```
세 번째 적분값:  5*Z/8
```
생각보다 깔끔하게 나왔네요? 이제 모든 항을 합치면 됩니다.

$$
\int \psi_t^* \hat{H} \psi_t \, d\tau = -Z^2 + (Z-2)(2Z) + \frac{5}{8}Z = Z^2 - \frac{27}{8}Z
$$

최소값을 찾으려면 미분해서 0이 되는 점을 찾으면 되죠. 이차식이니 간단합니다. $Z=\frac{27}{16}=1.6875$네요. 혹시 10장에서 유효 핵 전하를 얘기했을 때 나왔던 값, 기억하시나요? 같은 껍질의 $1s$ 전자끼리는 0.30만큼 핵 전하를 가려서 헬륨 원자의 유효 핵 전하는 1.70이라고 했었습니다. 그리고 방금 변분법으로 구한 $Z$ 값이 1.6875였으니 거의 비슷한 값입니다. Slater는 실험 데이터를 정리해서 만든 경험적인 규칙이었는데 실제로 변분법의 결과와 거의 같은 값을 주니 이론적인 값을 잘 근사하도록 맞춰준 것입니다.

이 값을 대입해보면 -2.847656 Hartree, 약 -77.49 eV가 나옵니다. 실험값과 약 1.5 eV 정도 차이 나네요. 오차가 작은 편은 아니지만 섭동론에서 1차 보정으로 계산한 값(-74.83 eV)보다는 훨씬 참값에 가깝습니다. 물론 섭동론은 2차 보정을 하면 변분법보다 더욱 정확도가 높은 값을 얻을 수 있죠. 그렇다고 해서 섭동론이 변분법보다 무조건 나은 방법이라고 할 수는 없습니다. 두 방법이 각자의 장단점이 있기 때문이죠.  


## 두 방법의 비교

변분법의 강점은 에너지의 상한이 보장된다는 것과 비교적 쉬운 계산 난이도입니다. 매개변수 하나만 조절했는데도 1.5 eV까지 오차를 줄였고, 참값보다 더 큰 값이라는 것을 알고 있으니 더 작은 값이 나온다면 더 정확한 값이라는 것을 바로 알 수 있습니다. 그리고 섭동론은 2차 보정부터는 손으로 계산할 수도 없는 수준으로 중간 과정이 복잡해지고요.
반면에 섭동론의 강점은 개선을 체계적으로 할 수 있다는 점입니다. 차수를 올릴수록 정확도가 개선되고 어떤 상태들이 섞이는지가 식에 드러나죠. 다만 참값보다 큰지 작은지는 알 수 없습니다.  

두 방법의 결과가 어떻게 다른지 그림으로도 살펴보겠습니다.
```python
import numpy as np
import matplotlib.pyplot as plt

Eh = 27.2114
E = lambda Z: (Z**2 - 27 / 8 * Z) * Eh
Zs = np.linspace(1.0, 2.6, 450)

fig, ax = plt.subplots(figsize=(7, 4))
ax.plot(Zs, E(Zs), lw=2, color="steelblue")
ax.plot(27 / 16, E(27 / 16), "o", ms=10, color="darkorange", label=f"variational")
ax.plot(2.0, E(2.0), "s", ms=9, color="crimson", label="Z=2 (perturbation 1st)")
ax.axhline(-79.00, color="black", ls="--", lw=1.5, label="experiment")
ax.set_xlabel("Z")
ax.set_ylabel("E (eV)")
ax.legend(fontsize=9)
ax.grid(alpha=0.3)
plt.show()
```
![변분법과 섭동론의 비교](/assets/image-88.png)

빨간색 점이 섭동론으로 구한 1차 보정 에너지인데 $Z=2$를 고정하고 에너지만 보정한 것이라 최소값과 좀 떨어져 있죠. 반면에 주황색 점은 변분법으로 구한 값으로 곡선의 최소점에 가 있습니다. 물론 그림에서는 섭동론의 2차 보정 항이 표현되어 있지 않지만 두 방법의 접근 방식이 어떻게 다른지는 알 수 있습니다.

그래서 실제로는 위에 팁 박스에서도 언급한 것처럼 둘을 섞어서 씁니다. 17장의 Hartree-Fock은 변분법인데 21장에서 볼 MP2는 그 위에 얹는 2차 섭동론이라 둘이 함께 쓰면 정확도가 좋아지거든요.


## 다음 이야기

변분법과 섭동론을 헬륨 원자에 적용해서 에너지를 얼마나 정확하게 얻을 수 있는지 알아보았습니다. 지금은 바닥 상태의 에너지를 계산하는 데에만 썼기 때문에 축퇴가 없는 경우여서 섭동론이 잘 맞았습니다. 하지만 11.3절에서 말했던 것처럼 축퇴가 있는 경우에는 지금과 같은 보정 식을 그대로 쓸 수는 없는데요, 다음 절에서 축퇴가 있는 경우에는 어떻게 처리하는지 살펴보도록 하겠습니다.


## 확인 문제
1. 헬륨 이온의 에너지는 -2 Hartree입니다. 위 결과를 이용해 헬륨의 첫 번째 이온화 에너지를 구하고 실제값 24.59 eV와 비교해보세요.
2. 리튬 이온(전자 두 개, $Z=3$)에 변분법을 적용해보세요. $E(Z) = Z^2 - \frac{43}{8}Z$가 나옵니다. 최적 $Z$값은 얼마인가요?
3. 변분법에서 얻은 유효 핵 전하 값 1.6875는 두 전자에 같은 값을 쓴 결과입니다. 두 전자에 서로 다른 유효 핵 전하 값을 써서 변분법을 적용하면 어떻게 될까요?
