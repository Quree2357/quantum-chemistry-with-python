# 1.1. 흑체 복사와 자외선 파국

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/01-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

## 흑체 복사

**흑체(blackbody)**란 들어오는 빛(전자기파)을 모두 흡수하는 가상의 물체로, 모든 빛을 흡수하여 검게 보이기 때문에 흑체(黑體)라는 이름이 붙었습니다. 이 흑체는 특이한 성질을 가지는데요, 바로 온도에 따라 특정한 파장의 빛을 방출한다는 것입니다. 주변에서 쉽게 볼 수 있는 예시로는 달군 쇠가 있겠네요. 쇠를 달구면 온도가 점점 올라가면서 빨간색, 노란색, 흰색과 같이 색을 내죠.

![달군 쇠](/assets/image-1.png)  
[출처: https://javalab.org/ko/stellar_classification/]

사실 별의 색깔도 흑체 복사와 같은 원리입니다. 별의 온도가 다르면 색이 다른 색으로 보이죠. 고등학교 때 배운 기억이 나실 겁니다.

![별의 온도와 색](/assets/image.png)  
[출처: https://astro.kasi.re.kr/learning/pageView/5042]

이렇게 흑체가 빛을 방출하는 현상을 **흑체 복사(blackbody radiation)**라고 합니다. 그런데 이게 왜 고전물리학에서 문제가 되었을까요?

## 고전물리학의 시도

당연히, 물리학자들은 이 현상을 고전물리학으로 설명하려고 노력했습니다. 레일리(Rayleigh)와 진즈(Jeans)라는 물리학자는 당시의 물리학 이론을 바탕으로 하여 흑체가 방출하는 빛의 스펙트럼을 계산하였습니다. 아래의 식이 바로 레일리-진즈 법칙으로 알려져 있는 식입니다.

$$
\rho_\nu(T)d\nu = \frac{8 \pi k_B T}{c^3} \nu^2 d\nu
$$

식이 조금 복잡하게 보일 수도 있겠지만, 천천히 따라가봅시다. 등호 왼쪽에서 $T$는 온도이고, $\nu$는 주파수(frequency)를 나타냅니다. 그러니까, 특정 온도에서 방출되는 특정 주파수의 에너지라고 말할 수 있겠네요. 오른쪽에 있는 $k_B$라는 값은 열역학 이론에서 사용되는 볼츠만 상수(Boltzmann constant)라고 불리는 상수입니다. 분모의 $c$는 빛의 속도를 뜻하는 상수입니다. 그러니까 특정 온도를 고정해놓으면 분수로 나타낸 값은 모두 상수라고 봐도 되겠네요. 그리고 그 오른쪽에 주파수의 제곱 항이 있습니다. 다시 말하면, 흑체가 방출하는 빛의 에너지는 주파수가 증가함에 따라 제곱으로 증가한다는 뜻입니다. 한번 그림으로 그려볼까요?

```python
import numpy as np
import matplotlib.pyplot as plt
from scipy.constants import k, c

T = 6000
nu = np.linspace(1e12, 1e15, 2000)


def rayleigh_jeans(nu, T):
    return 8 * np.pi * k * T / c**3 * nu**2


plt.plot(nu, rayleigh_jeans(nu, T))
plt.xlabel(r"Frequency  $\nu$  [Hz]")
plt.ylabel(r"Energy density")
plt.xlim(0, 1e15)
plt.ylim(0, 8e-14)
plt.title(f"Blackbody radiation at T = {T} K")
plt.show()
```
![레일리-진즈 법칙](/assets/image-2.png)  

이 모델의 문제가 무엇인지 아시겠나요? 실제 실험 데이터를 같이 표시하지는 않았지만 주파수가 커질수록 에너지가 급격히 커지고 있습니다. 이 식대로라면 우리 우주를 이루는 물질은 모두 순식간에 모든 에너지를 자외선과 X선으로 방출하고 얼어붙어야 합니다! 이를 **자외선 파국(ultraviolet catastrophe)**이라고 부릅니다. 하지만 실제 우리 우주는 그렇지 않습니다. 분명히 레일리와 진즈는 고전 이론에 충실하게 공식을 유도하였고, 중간에 틀린 점도 없었습니다. 그렇다면 대체 무엇이 문제였을까요?

## 플랑크의 양자 가설

여기서 플랑크(Planck)라는 물리학자가 등장합니다. 이름은 다들 익숙하실 텐데요. 1900년, 플랑크는 흑체 복사를 설명하기 위해서 아주 특이한 가정을 하나 합니다. 바로, 에너지가 연속적이지 않다는 것입니다. 이게 대체 무슨 소리지 싶으시겠지만 자세히 설명드리겠습니다.  

플랑크는 흑체 복사를 흑체의 벽면을 이루는 입자들이 진동하면서 빛을 내는 것이라고 생각했습니다. 그러니까 흑체가 **진동자(oscillator)**들로 이루어져 있다는 말이죠. 그러니까 이 진동자들의 진동수에 따라 서로 다른 주파수의 빛이 나온다는 것이죠. 예를 들면 1초에 100번 진동하는 진동자는 100Hz의 빛을 내고, 1초에 200번 진동하는 진동자는 200Hz의 빛을 내고, 이런 식입니다. 여기에서 플랑크는 아까 말한 특이한 가정을 추가하였습니다. 진동자가 가질 수 있는 에너지는 아무 값이나 가능한 것이 아니라, 정해진 크기의 덩어리 단위로만 에너지를 가질 수 있다는 것입니다. 정확히는 $E=nh\nu$라는 식으로 표현했습니다. 여기서 $E$는 전자의 에너지를 나타내고, $\nu$는 아까처럼 주파수를 나타냅니다. $n$은 정수인데요, 그러니까 0, 1, 2, 3과 같은 수입니다. 다시 말하면, 진동자가 가질 수 있는 에너지는 $h\nu$의 정수배밖에 가질 수 없다는 뜻입니다. 그리고 이 식에서 $h$가 바로 플랑크 상수라고 부르는 녀석입니다.

...왜요? 라고 묻기 전에 일단 이 가정으로부터 출발하여 유도된 플랑크의 흑체 복사 모델을 살펴봅시다. 플랑크는 다음과 같은 식을 유도하였습니다.  

$$
\rho_\nu(T)d\nu = \frac{8 \pi h}{c^3} \frac{\nu^3}{e^{h\nu/k_BT}-1} d\nu
$$

아까의 레일리-진즈 식보다 훨씬 더 복잡해졌네요! 여기서 $h$는 플랑크 상수라고 합니다. 일단 그림을 먼저 그려봅시다.

```python
from scipy.constants import h

T = 6000
nu = np.linspace(1e12, 1e15, 2000)


def planck(nu, T):
    return 8 * np.pi * h / c**3 * nu**3 / np.expm1(h * nu / (k * T))


plt.plot(nu, planck(nu, T))
plt.xlabel(r"Frequency  $\nu$  [Hz]")
plt.ylabel(r"Energy density")
plt.xlim(0, 1e15)
plt.ylim(0, 2e-15)
plt.title(f"Blackbody radiation at T = {T} K")
plt.show()
```

![플랑크의 법칙](/assets/image-3.png)  

우와, 아까랑 전혀 다른 모양의 그래프가 나왔습니다! 더구나 이제는 진동수가 커진다고 해서 에너지가 무한대로 발산하지도 않습니다. 그리고 더 놀라운 사실은, 이 그래프는 실험으로 얻은 데이터와 정확히 일치했습니다. 플랑크가 흑체 복사 스펙트럼을 설명하는 데 성공한 것입니다! 이상한 가정을 하나 추가하였을 뿐인데 모델이 전혀 달라졌습니다. 이 가정을 **양자 가설(quantum hypothesis)**이라고 부릅니다. 앞의 레일리-진즈 식과 비교해볼까요?  
```python
plt.plot(nu, rayleigh_jeans(nu, T), label="Raylegh-Jeans")
plt.plot(nu, planck(nu, T), label="Planck")
plt.xlabel(r"Frequency  $\nu$  [Hz]")
plt.ylabel(r"Energy density")
plt.xlim(0, 1e15)
plt.ylim(0, np.max(planck(nu, T)) * 2)
plt.title(f"Blackbody radiation at T = {T} K")
plt.legend()
plt.show()
```
![레일리-진즈 식과 플랑크 식의 비교](/assets/image-19.png)  
한 가지 주목할 점은 두 식 모두 매우 낮은 주파수에서는 정확히 일치하고 있다는 점입니다. 그래서 당시 물리학자들이 더욱 진지하게 받아들였고요.

온도가 달라지면 어떻게 될까요?
```python
T = [3000, 4000, 5000, 6000]

for temperature in T:
    plt.plot(nu, planck(nu, temperature), label=f"T = {temperature} K")
plt.xlabel(r"Frequency  $\nu$  [Hz]")
plt.ylabel(r"Energy density")
plt.xlim(0, 1e15)
plt.ylim(0, 2e-15)
plt.title("Blackbody radiation at different temperatures")
plt.legend()
plt.show()
```
![다양한 온도에서의 흑체 복사](/assets/image-4.png)  

온도가 높아질수록 흑체 복사의 에너지가 증가합니다. 또 하나 중요한 점은 에너지가 최대인 지점의 진동수 또한 증가하고 있다는 것입니다. 진동수가 증가한다는 것은 반대로 파장이 감소하고 있다는 뜻인데요, 온도가 증가함에 따라 긴 파장인 빨간색에서 짧은 파장인 파란색으로 가고 있습니다. 처음에 얘기한 별의 온도에 따른 색깔과 같은 결론이네요! 역시 플랑크의 흑체 복사 모델이 맞는 것 같습니다.

하지만 당시의 물리학자들은 플랑크가 추가한 가정을 받아들이지 않았습니다. 에너지가 연속적이지 않다뇨! 그럴 리가 없었습니다. 무언가 잘못된 것이라고 생각했습니다. 그러나 몇 년 뒤, 상황은 바뀌었습니다.
