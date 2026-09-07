# 11.1. 변분 원리

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/11-01.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

먼저 변분법에 대해서 살펴봅시다. Euler와 Lagrange에 의해 정립된 변분법은 양자역학뿐만 아니라 수학과 물리학 전반에서 쓰이는 방법인데, 함수의 최대/최소값을 찾는 방식입니다. 5.3절에서 이런 이야기를 잠깐 하고 넘어갔었죠.

> 아무 함수나 가져와 에너지 기댓값을 계산해도 이 부등식이 성립합니다. (...) 이 성질을 변분 원리(variational principle)라고 부르는데, (...) 11장에서 제대로 다루겠습니다.

변분법은 정답을 몰라도 근사의 정확도를 채점할 수 있다는 점에서 매우 유용한 방법입니다. 한번 자세히 들여다보죠.


##
