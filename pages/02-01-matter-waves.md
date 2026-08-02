# 2.1. 물질도 파동이다

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/02-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

## 질문을 뒤집다

1924년, 프랑스의 de Broglie는 박사 학위 논문을 제출합니다. 분량이 얇았고, 심사위원들은 이걸 어떻게 판단해야 할지 몰라 곤란해했다고 합니다. 결국 Einstein에게 의견을 물었고, Einstein은 "그는 거대한 장막의 한 귀퉁이를 들어 올렸다"고 답했습니다.

de Broglie가 한 일은 질문을 거꾸로 생각해본 것입니다. 빛을 생각해봅시다. 200년 넘게 파동이라고 믿어왔는데, 광전 효과에서 입자처럼 행동했습니다. 그러니까 빛은 파동인 줄 알았는데 입자성도 가지고 있었던 겁니다. 그러면 반대로, 입자인 줄 알았던 전자도 사실은 파동성을 가지고 있는 것 아닐까요? 이상하게 들릴 수도 있겠지만 그런 생각을 가질 만합니다.


## de Broglie의 물질파

전자가 파동이라면 파장이 있어야겠죠. 그런데 대체 얼마일까요? de Broglie는 빛에서 힌트를 얻었습니다. 광자의 에너지는 $E = h\nu$이고, 상대성 이론에 따르면 질량이 없는 입자의 에너지와 운동량은 $E = pc$의 관계에 있습니다. 두 식을 합치면 다음과 같이 되겠네요.

$$
p = \frac{h\nu}{c} = \frac{h}{\lambda}
$$

즉, 광자의 파장과 운동량은 $\lambda = h/p$라는 식으로 이어져 있습니다. 여기서 de Broglie는 대담한 도약을 합니다. 이 관계가 빛에만 적용될 이유가 있을까? 입자도 파동성을 나타낸다면 이 식을 따르지 않을까? 속도 $v$로 운동하는 질량이 $m$인 입자의 운동량이 $p=mv$로 주어진다는 것은 잘 알고 계실 겁니다.

$$
\lambda = \frac{h}{p} = \frac{h}{mv}
$$

이것이 바로 **de Broglie 파장(de Broglie wavelength)**입니다. 질량과 속도를 가진 모든 물체에 파장이 있다는 주장이죠.

식만 놓고 보면 별것 아닌 것 같지만, 주장하는 바는 어마어마합니다. 여러분이 지금 움직이고 있다면 여러분에게도 파장이 있다는 뜻이니까요.


## 그런데 왜 아무도 눈치채지 못했을까

당연히 이런 반론이 나옵니다. 물질이 파동이라면 왜 우리는 한 번도 그걸 본 적이 없을까요? 야구공이 회절하는 걸 본 사람은 아무도 없잖아요. 계산해보면 답이 나옵니다.

```python
import numpy as np
from scipy.constants import h


def de_broglie(m, v):
    return h / (m * v)


for name, m, v in [("야구공", 0.1, 50), ("걷는 사람", 70, 1), ("자동차", 1000, 30)]:
    print(f"{name:10s}  lambda = {de_broglie(m, v):.3e} m")
```
```
=== 일상 물체 ===
야구공         lambda = 1.325e-34 m
걷는 사람       lambda = 9.466e-36 m
자동차         lambda = 2.209e-38 m
```

야구공의 파장은 $10^{-34}$ m입니다. 원자 크기가 $10^{-10}$ m, 원자핵이 $10^{-15}$ m인데, 야구공의 파장은 원자핵의 크기보다도 자릿수가 19개나 작습니다.

파동의 성질이 드러나려면 파장과 비슷한 크기를 가진 구조를 통과해야 합니다. 소리가 문틈을 돌아 나오는 것도 소리의 파장이 문틈 크기와 비슷하기 때문이죠. 그런데 $10^{-34}$ m짜리 틈은 이 우주에서 만들 수 없습니다. 야구공은 파동성을 가지지 않는 것이 아니라 그 성질을 관측할 방법이 없는 것입니다.

그러면 전자는 어떨까요? 전자의 질량은 $9.1\times10^{-31}$ kg입니다. 대충 광속의 1%의 속도로 움직이고 있다고 가정해보죠.

```python
from scipy.constants import c, m_e

print(f"3000km/s 속도로 움직이는 전자의 파장 = {de_broglie(m_e, c * 0.01):.3e} m")
```
```
3000km/s 속도로 움직이는 전자의 파장 = 2.426e-10 m
```

전자의 파장은 약 0.24 nm입니다. 원자의 크기 정도네요. 원자가 규칙적으로 배열된 것이 바로 결정이죠. 결정 속 원자 사이 간격이 딱 이 정도입니다. 그러니까 전자를 결정에 쏘면 회절이 일어나야 합니다. X선 회절과 똑같은 원리로요.

그림으로 정리하면 이렇습니다.

```python
from scipy.constants import e
import matplotlib.pyplot as plt

items = [
    ("electron", m_e * (c * 0.01), "red"),
    ("He atom (300 K)", 6.6e-27 * 1.25e3, "orange"),
    ("baseball", 0.1 * 50, "green"),
    ("human (walking)", 70 * 1, "blue"),
]

fig, ax = plt.subplots(figsize=(8, 4.5))
for i, (name, p, col) in enumerate(items):
    lam = h / p
    ax.barh(i, lam, color=col, height=0.55)
    ax.text(lam * 1.6, i, f"{lam:.1e} m", va="center", fontsize=9)

ax.axvline(1e-10, color="black", ls="--", lw=1.5)
ax.axvline(1e-15, color="gray", ls=":", lw=1.5)

ax.set_xscale("log")
ax.set_yticks(range(len(items)))
ax.set_yticklabels([n for n, _, _ in items], fontsize=9)
ax.set_xlim(1e-38, 1e-6)
ax.set_xlabel("de Broglie wavelength (m)")
ax.grid(alpha=0.3, axis="x", which="both")
plt.show()
```
![물질파 파장 비교](/assets/image-26.png)

굵은 점선이 원자의 크기이고 얇은 점선이 원자핵의 크기입니다. 전자와 헬륨 원자 정도만 그 근처에 있고 나머지는 까마득히 왼쪽에 있죠. de Broglie의 주장을 검증하려면 전자를 써야 한다는 것이 이 그림의 결론입니다.


## 전자도 회절한다

1924년, 미국의 Davisson과 Germer가 니켈 결정에 전자를 쏘는 실험을 하고 있었습니다. 원래는 니켈의 표면을 연구하는 목적의 실험이었는데, 사고로 진공관이 깨져 니켈이 산화되었고, 이를 되살리려고 고온으로 가열하는 과정에서 니켈이 큰 단결정으로 자라버렸습니다. 그리고 그 결정에 전자를 쏘자 특정 각도에서만 전자가 유난히 많이 튀어나왔습니다. 54 V로 가속한 전자(즉, 운동 에너지가 54eV인 전자)를 쏘았을 때 50도 각도 방향에서 뚜렷한 봉우리가 있는 회절 무늬가 나타났습니다.
Davisson은 연구 도중에 de Broglie의 공식을 알게 되었고 이 공식으로 실험 결과를 설명하려고 했습니다. 그리고 Davisson과 Germer는 1927년 발표한 논문에서 니켈 결정 실험으로 얻은 데이터와 de Broglie 공식이 일치한다는 것을 보였습니다.

정말 de Broglie 공식과 맞는지 확인해봅시다. 파장을 두 가지 방법으로 구해서 비교하면 됩니다. 회절 무늬로부터 파장을 구하는 방법은 Bragg 법칙를 이용하면 됩니다.

```python
# 방법 1: de Broglie 관계식으로 계산
def de_broglie_electron(V):
    """전위차 V로 가속한 전자의 de Broglie 파장"""
    p = np.sqrt(2 * m_e * e * V)  # 전자의 운동에너지 (eV) = p^2/2m
    return h / p


lam_deBroglie = de_broglie_electron(54)

# 방법 2: 회절 무늬에서 읽기
d = 215e-12  # 니켈의 원자면 간격 (m)
theta = np.deg2rad(50)  # 봉우리가 나타난 각도
lam_diffraction = d * np.sin(theta)

print(f"de Broglie 관계식  : {lam_deBroglie * 1e12:.1f} pm")
print(f"회절 무늬에서    : {lam_diffraction * 1e12:.1f} pm")
print(f"오차            : {abs(lam_deBroglie - lam_diffraction) / lam_diffraction * 100:.1f} %")
```
```
de Broglie 관계식  : 166.9 pm
회절 무늬에서    : 164.7 pm
차이            : 1.3 %
```

거의 일치하네요! 두 값이 완전히 다른 경로에서 나왔다는 점이 중요합니다. 왼쪽은 Planck 상수와 가속 전압만으로 계산한 값이고, 오른쪽은 니켈의 결정 구조와 관측된 각도에서 읽은 값입니다. 서로 아무 관계도 없어 보이는 두 계산이 같은 숫자를 내놓았죠.

같은 해 영국에서는 G. P. Thomson이 얇은 금속박에 전자를 통과시켜 동심원 회절 무늬를 얻었습니다. 여기에 재미있는 사실이 하나 있습니다. G. P. Thomson의 아버지가 바로 1897년에 전자를 발견한 J. J. Thomson입니다. 아버지는 전자가 입자임을 밝혀 1906년에 노벨상을 받았고, 아들은 전자가 파동임을 밝혀 1937년에 노벨상을 받았습니다. 

de Broglie는 1929년 노벨 물리학상을 받습니다. 박사 논문을 낸 지 5년 만이었습니다. 


## 정리하며

빛이 파동이면서 입자성을 갖는다면, 물질도 입자이면서 파동성을 가질 수 있다는 것이 de Broglie의 발상이었습니다. de Broglie 파장은 $\lambda = h/p$입니다. 광자에 성립하던 관계를 물질로 확장한 것이죠. 일상 물체의 파장은 원자핵보다도 작아서 관측할 방법이 없습니다. 파동성이 없는 것이 아니라 드러날 기회가 없는 것입니다.

그런데 이 아이디어에는 실험 확인보다 더 놀라운 결과가 딸려 있었습니다. 1장에서 Bohr가 설명하지 못하고 남겨둔 그 가정, **각운동량이 왜 $\hbar$의 정수배인가**에 대한 답이 여기서 저절로 나옵니다. 다음 절에서 이야기를 계속해보죠.


## 확인 문제

1. 전자를 10000 V로 가속하면 파장이 얼마가 될까요? 이 정도 파장이면 무엇을 관찰할 수 있을까요? (힌트: 전자현미경)
2. 중성자의 질량은 전자의 약 1840배입니다. 같은 속도로 움직이는 중성자와 전자의 de Broglie 파장은 몇 배 차이가 날까요?
3. 파장이 0.1 nm인 전자를 얻으려면 몇 볼트로 가속해야 할까요? 위의 `de_broglie_electron` 함수를 거꾸로 풀어보세요.