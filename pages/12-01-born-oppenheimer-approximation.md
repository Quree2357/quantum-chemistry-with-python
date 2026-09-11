# 12.1. Born-Oppenheimer 근사



드디어 분자를 다뤄볼 차례입니다. 언제나 시작은 Schrödinger 방정식을 만드는 것입니다. 상황을 단순하게 만들기 위해 원자핵과 전자는 모두 크기가 없는 점이고, 스핀-궤도 결합이나 상대론적인 상호작용 같은 건 모두 무시한다고 가정해봅시다. 그러면 분자의 Hamiltonian은 (원자 단위계로) 이렇게 쓸 수 있습니다.

$$
\hat{H} = \underbrace{-\frac{1}{2}\sum_{a} \frac{1}{m_a} \nabla_a^2}_{\text{a번 핵의 운동에너지}} \overbrace{-\frac{1}{2}\sum_{i} \frac{1}{m_e} \nabla_i^2}^{\text{i번 전자의 운동에너지}} + \underbrace{\sum_a \sum_{b>a}\frac{Z_a Z_b}{r_{ab}} - \sum_a \sum_{i}\frac{Z_a}{r_{ai}} + \sum_i \sum_{j>i}\frac{1}{r_{ij}}}_{\text{핵/핵, 핵/전자, 전자/전자간 반발 항}}
$$

항이 무지막지하게 늘었습니다. 입자가 3개만 있어도 정확한 해석해를 못 구한다고 했는데 큰일났네요. 이 말도 안 되는 몹을 어떻게 잡아야 한담...


## 너무 무거워서 안 움직여요

여기 좋은 방법이 하나 있습니다. 원자핵에 들어있는 양성자나 중성자는 전자보다 약 1840배 정도 무겁습니다. 9장에서 환산질량을 계산할 때는 이 사실을 전자의 질량으로 근사하는 데 썼었죠. 이번에는 조금 다른 관점에서 봅시다.  
원자핵과 전자는 정전기적 인력으로 서로를 붙잡고 공간 속을 떠다니고 있을 겁니다. 그렇다면 외부에서 알짜힘이 가해지지 않는 한 질량중심은 그대로 있겠죠. 그런데 원자핵이 전자보다 엄청나게 무거우니 사실상 원자핵은 가만히 있고 전자만 움직인다고 봐도 될 겁니다.

이렇게 보는 것을 **Born-Oppenheimer 근사**라고 합니다. 2.3절에서 잠시 만났던 Born과 그의 대학원생이었던 Oppenheimer가 1927년에 제안한 방법이죠.(영화로 제작된 그 사람이 맞습니다.)  
핵들의 위치를 고정시키고 나면 핵의 운동 에너지에 해당하는 항이 사라지고, 핵/핵 반발 항은 그냥 상수 $V_{NN}$이 되어버립니다.

$$
\hat{H_e} = -\frac{1}{2}\sum_{i} \frac{1}{m_e} \nabla_i^2 + V_{NN} - \sum_a \sum_{i}\frac{Z_a}{r_{ai}} + \sum_i \sum_{j>i}\frac{1}{r_{ij}}
$$

이제 전자만의 문제가 되었습니다. 이 식을 풀면 전자의 파동 함수와 에너지가 나오겠죠. 물론, 핵들의 가능한 배치는 무한히 많고, 각각의 배치마다 서로 다른 식과 해가 나올 겁니다. 그래서 전자의 에너지와 파동 함수는 핵들의 위치 $\mathbf{R}$에도 의존합니다.

$$
\hat{H_e}(\mathbf{R})\psi_e(r;\mathbf{R}) = E_e(\mathbf{R})\psi_e(r;\mathbf{R})
$$

전체 에너지를 구하려면 $E_e(\mathbf{R})$에 핵/핵 반발 항을 더해주면 되죠. 전체 파동 함수를 구하려면 $\psi_e(r;\mathbf{R})$에 핵의 파동 함수를 곱해주면 되고요.

$$
E(\mathbf{R}) = E_e(\mathbf{R}) + V_{NN}
$$
$$
\Psi(r;\mathbf{R})=\psi_e(r;\mathbf{R})\,\chi_N(\mathbf{R})
$$


## 이제 내가 움직일 차례인가?

이제 여기서 재미있는 일이 벌어집니다. 위의 식을 잘 보면 전체 에너지가 핵들의 위치 $\mathbf{R}$에만 의존하죠. 그러니까 에너지를 $\mathbf{R}$에 대한 함수로 생각할 수 있습니다.  
이 함수 $E(\mathbf{R})$을 **퍼텐셜 에너지 곡면(potential energy surface, PES)**이라고 부르는데요. 전자의 에너지가 핵 입장에서 보면 퍼텐셜로 느껴지기 때문에 이런 이름이 붙었습니다. 결국 퍼텐셜 에너지 곡면에서 최소점에 해당하는 핵의 위치가 가장 안정한 구조라는 뜻이죠.  

핵들의 위치에 따라 분자의 전체 에너지가 어떻게 변하는지를 나타내는 것이 PES이니 이 곡면을 따라서 분자의 구조 완화, 혹은 분자 간 화학 반응의 경로가 정해지게 됩니다. 계산화학 프로그램에서 분자의 기하 최적화 같은 계산을 돌릴 때 이런 과정을 거치는 거죠.

이산화탄소처럼 3개의 원자로 이루어진 선형 분자의 경우에는 PES를 그리면 대충 이런 식으로 그려지게 됩니다.
```python
import numpy as np
import matplotlib.pyplot as plt

De, a, r0 = 11.372, 2.075, 1.173
r1 = np.linspace(0.8, 4, 150)
r2 = np.linspace(0.8, 4, 150)

X, Y = np.meshgrid(r1, r2)
V = De * ((1 - np.exp(-a * (X - r0))) ** 2 + (1 - np.exp(-a * (Y - r0))) ** 2)

fig = plt.figure(figsize=(12, 6))

ax1 = fig.add_subplot(1, 2, 1, projection="3d")
ax1.plot_surface(X, Y, V, cmap="viridis", rstride=3, cstride=3, linewidth=0, antialiased=True, alpha=0.9)
ax1.contour(X, Y, V, levels=12, zdir="z", offset=V.min() - 2, cmap="viridis", linewidths=1)
ax1.plot([r0], [r0], [V.min()], "o", ms=8, color="crimson")
ax1.set_xlabel("R1 (Angstrom)")
ax1.set_ylabel("R2 (Angstrom)")
ax1.set_zlabel("E (eV)")
ax1.set_xlim(0.8, 4)
ax1.set_ylim(0.8, 4)
ax1.set_zlim(-2, 22)
ax1.view_init(elev=28, azim=-110)

ax2 = fig.add_subplot(1, 2, 2)
cs = ax2.contourf(X, Y, V, levels=25, cmap="viridis")
ax2.contour(X, Y, V, levels=12, colors="white", linewidths=1, alpha=0.5)
ax2.plot(r0, r0, "o", ms=10, color="crimson")
ax2.set_xlabel("R1 (Angstrom)")
ax2.set_ylabel("R2 (Angstrom)")
plt.colorbar(cs, ax=ax2, label="E (eV)")

plt.tight_layout()
plt.show()
```
![이산화탄소의 단순 PES](/assets/image-89.png)

탄소와 산소 사이의 결합 길이가 달라지면 핵들의 상대적인 위치가 달라지면서 분자의 전체 에너지가 바뀌죠. 빨간 점으로 표시한 부분이 가장 에너지가 낮아 가장 안정한 구조의 이산화탄소 분자입니다. 이 곡면에 공을 굴린다고 생각해보세요. 공은 결국 가장 낮은 곳으로 내려갈 거고, 그 지점이 가장 안정한 곳인거죠.

$N$개의 원자로 이루어진 분자라면 각각의 원자핵이 3차원 공간 상에서 자유롭게 움직일 수 있으니 PES를 그리려면 $3N$개의 변수가 필요할 겁니다. (분자 전체가 통째로 병진 운동이나 회전 운동을 하는 경우를 제외하면 전체 자유도(변수의 개수)는 $3N-6$ 입니다.) 아쉽게도 그림으로 그려볼 수 있는 건 여기까지겠네요...
