# 14.4. 물 분자의 분자 오비탈

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/14-04.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

세 개 이상의 원자로 이루어진 분자도 마찬가지 방법으로 분자 오비탈을 구성할 수 있습니다. 원자 오비탈들을 에너지와 대칭을 고려해서 잘 섞으면 되죠. 어떤 오비탈을 어떻게 섞어야 하는지는 Hamiltonian 행렬이 어떻게 쪼개지는지를 보면 알 수 있습니다. 물 분자를 예로 들어보죠.


## 물 분자 오비탈의 섞임

13.4절에서 만들었던 물 분자의 Hamiltonian 행렬을 다시 봅시다. 그때는 산소 원자의 $1s$ 오비탈까지 고려했었지만 14.2절에서 얘기했던 것처럼 에너지 차이를 생각하면 자기들끼리만 섞일 겁니다. 산소의 $1s$ 오비탈을 제외한 나머지 여섯 개만 보도록 하죠.
```
              O 2s  O 2pz     H+  O 2px  O 2py     H-
  O 2s (A1)      x      x      x      .      .      .
 O 2pz (A1)      x      x      x      .      .      .
    H+ (A1)      x      x      x      .      .      .
 O 2px (B1)      .      .      .      x      .      .
 O 2py (B2)      .      .      .      .      x      x
    H- (B2)      .      .      .      .      x      x
```
참고로, 수소의 SALC는 13.5절에서 만들었던 것과 같은데 여기서는 $S=0$으로 두었기 때문에 규격화 상수가 $\frac{1}{2}$가 아니라 $\frac{1}{\sqrt{2}}$가 됩니다.  

산소의 $2p_x$ 오비탈은 13.4절에서 얘기했던 것처럼 대칭 종 때문에 섞이지 않고 비결합성 오비탈로 남습니다. 나머지 블록에서는 수소와 산소의 오비탈들이 섞이죠. 이때, 주의해야 할 것이 있습니다. 우리는 물 분자가 어떻게 생겼는지 이미 알고 있지만 원자들로부터 만들어질 때는 아무것도 모르는 상태에서 결합한다는 것입니다. 즉, 물 분자의 결합각이 104.5˚라는 사실을 모른다고 가정하는 겁니다.  

만약 결합각을 $\theta$라고 하고 산소의 $2p$ 오비탈이 수소의 $1s$ 오비탈과 정면으로 겹칠 때의 공명 적분을 $\beta_{sp}$라고 하면, 산소의 $2p_z$ 오비탈은 $\beta_{sp} \cos(\theta/2)$만큼만, $2p_y$ 오비탈은 $\beta_{sp} \sin(\theta/2)$만큼만 상호작용할 겁니다. 이번에도 $\beta$가 음수가 되도록 $2p_z$ 오비탈은 수소 쪽($-z$ 방향)을, $2p_y$ 오비탈은 $H_1$ 수소 방향이 되도록 잡겠습니다. 그러면 각 공명 적분은 이렇게 쓸 수 있습니다.

$$
\left\langle 2s_{\text{O}} \middle| \hat{H} \middle| \text{H}_{+} \right\rangle = \sqrt{2} \,\beta_{sp}, \quad
\left\langle 2p_{z,\text{O}} \middle| \hat{H} \middle| \text{H}_{+} \right\rangle = \sqrt{2} \,\beta_{sp} \cos \left( \frac{\theta}{2} \right), \quad
\left\langle 2p_{y,\text{O}} \middle| \hat{H} \middle| \text{H}_{-} \right\rangle = \sqrt{2} \,\beta_{sp} \sin \left( \frac{\theta}{2} \right)
$$


## 분자 오비탈 만들기

이제 세 블록에 대해 풀면 됩니다. $\alpha$는 14.2절에서 썼던 원자가 오비탈 에너지를, $\beta$는 대략적인 값을 쓰겠습니다. 결합각에 따라 값이 달라질테니 결합각 값을 입력받을 수 있도록 함수를 짜보죠.

```python
import numpy as np

a_O2s, a_O2p, a_H = -32.38, -15.85, -13.61
b_s, b_p = -4.0, -5.0


def water_blocks(angle):
    theta = np.radians(angle) / 2
    H_A1 = np.array(
        [
            [a_O2s, 0.0, np.sqrt(2) * b_s],
            [0.0, a_O2p, np.sqrt(2) * b_p * np.cos(theta)],
            [np.sqrt(2) * b_s, np.sqrt(2) * b_p * np.cos(theta), a_H],
        ]
    )
    H_B2 = np.array(
        [
            [a_O2p, np.sqrt(2) * b_p * np.sin(theta)],
            [np.sqrt(2) * b_p * np.sin(theta), a_H],
        ]
    )
    H_B1 = np.array([[a_O2p]])
    return {"a1": H_A1, "b2": H_B2, "b1": H_B1}


first = {"a1": 2, "b2": 1, "b1": 1}

levels = []
for sym, H in water_blocks(104.5).items():
    E, C = np.linalg.eigh(H)
    for i in range(len(E)):
        w_H = C[-1, i] ** 2 if sym != "b1" else 0.0
        levels.append((E[i], f"{first[sym] + i}{sym}", w_H))
levels.sort()

IE = {"2a1": 32.2, "1b2": 18.5, "3a1": 14.7, "1b1": 12.6}


n = 8
print(f"{'MO':>4} {'E (eV)':>8} {'H weight':>9} {'e':>3} {'IE (eV)':>9}")
for E, name, w_H in levels:
    k = min(n, 2)
    n -= k
    ie = f"{IE[name]:9.1f}" if name in IE else f"{'':>9}"
    print(f"{name:>4} {E:8.2f} {w_H:9.1%} {k:3d} {ie}")
```
```
  MO   E (eV)  H weight   e   IE (eV)
 2a1   -34.03      7.8%   2      32.2
 1b2   -20.43     40.2%   2      18.5
 3a1   -18.47     25.6%   2      14.7
 1b1   -15.85      0.0%   2      12.6
 4a1    -9.34     66.5%   0          
 2b2    -9.03     59.8%   0
```

일단 실제 결합각인 104.5˚에 대해 계산을 해봤습니다. 하나씩 살펴볼까요?  

$2a_1$ MO는 대부분이 산소의 $2s$ 오비탈로 이루어져 있네요. 수소의 $1s$ 오비탈과 에너지 차이가 약 20 eV나 되어 거의 섞이지 않았습니다. $1b_2$ 오비탈과 $3a_1$ 오비탈은 각각 산소의 $2p_y$와 $2p_z$ 오비탈이 섞인 것인데, $1b_2$의 에너지가 조금 더 낮고 수소 오비탈의 비중도 더 큽니다. 결합각 $\theta$가 90˚보다 크기 때문에 $\sin(\theta/2)$가 $\cos(\theta/2)$보다 큰 값을 가지기 때문이죠. 마지막으로 $1b_1$ 오비탈은 위에서 봤던 것처럼 산소의 $2p_x$ 오비탈이 그대로 비결합성으로 남은 겁니다.

그림으로 그리면 이렇게 되겠네요.

(그림 그리고 코드 지우기)
```python
import matplotlib.pyplot as plt

H2O_MO_LEVELS = [
    {"key": "2a1", "label": r"$2a_1$", "energy": 1.2, "symmetry": "a1", "type": "bonding"},
    {"key": "1b2", "label": r"$1b_2$", "energy": 3.9, "symmetry": "b2", "type": "bonding"},
    {"key": "3a1", "label": r"$3a_1$", "energy": 4.5, "symmetry": "a1", "type": "nonbonding"},
    {"key": "1b1", "label": r"$1b_1$", "energy": 4.8, "symmetry": "b1", "type": "nonbonding"},
    {"key": "4a1", "label": r"$4a_1$", "energy": 7.5, "symmetry": "a1", "type": "antibonding"},
    {"key": "2b2", "label": r"$2b_2$", "energy": 8.0, "symmetry": "b2", "type": "antibonding"},
]


def draw_energy_line(ax, x_center, y, width=0.8, color="black", linewidth=3, zorder=3):
    ax.plot(
        [x_center - width / 2, x_center + width / 2],
        [y, y],
        color=color,
        linewidth=linewidth,
        solid_capstyle="butt",
        zorder=zorder,
    )


def draw_orbital_group(
    ax, x_center, y, number_of_orbitals=1, orbital_spacing=0.60, level_width=0.46, level_color="black"
):

    if number_of_orbitals == 1:
        x_positions = [x_center]
    else:
        start = x_center - orbital_spacing * (number_of_orbitals - 1) / 2
        x_positions = [start + i * orbital_spacing for i in range(number_of_orbitals)]

    for x in x_positions:
        draw_energy_line(ax, x_center=x, y=y, width=level_width, color=level_color)

    return x_positions


def draw_correlation_line(ax, x1, y1, x2, y2, color="0.75", linewidth=1.0, linestyle="--"):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle=linestyle, zorder=1)


def connect_from_left(ax, atomic_x, atomic_y, molecular_x, molecular_y, color="0.75"):
    draw_correlation_line(ax, x1=atomic_x + 0.23, y1=atomic_y, x2=molecular_x - 0.31, y2=molecular_y, color=color)


def connect_from_right(ax, fragment_x, fragment_y, molecular_x, molecular_y, color="0.75"):
    draw_correlation_line(ax, x1=fragment_x - 0.23, y1=fragment_y, x2=molecular_x + 0.31, y2=molecular_y, color=color)


def get_mo_level(key):
    for level in H2O_MO_LEVELS:
        if level["key"] == key:
            return level

    raise KeyError(f"Unknown MO key: {key}")


def get_mo_energy(key):
    return get_mo_level(key)["energy"]


def get_level_color(level_type):
    colors = {"bonding": "tab:blue", "nonbonding": "black", "antibonding": "tab:red"}

    return colors[level_type]


def plot_h2o_mo_diagram(show_correlation_lines=True):
    x_o = 1.5
    x_mo = 5.5
    x_h2 = 10.0

    y_o_2s = 2.0
    y_o_2p = 4.8

    y_h_salc = 4.3

    fig, ax = plt.subplots(figsize=(9, 10))

    o_2s_positions = draw_orbital_group(ax, x_center=x_o, y=y_o_2s, number_of_orbitals=1, level_color="tab:green")
    o_2p_positions = draw_orbital_group(
        ax, x_center=x_o, y=y_o_2p, number_of_orbitals=3, orbital_spacing=0.7, level_color="tab:green"
    )

    o_2px = o_2p_positions[0]
    o_2py = o_2p_positions[1]
    o_2pz = o_2p_positions[2]

    ax.text(o_2px, y_o_2p + 0.25, r"$2p_x$", fontsize=12, ha="center")
    ax.text(o_2py, y_o_2p + 0.25, r"$2p_y$", fontsize=12, ha="center")
    ax.text(o_2pz, y_o_2p + 0.25, r"$2p_z$", fontsize=12, ha="center")
    ax.text(x_o - 0.65, y_o_2s, r"$2s$", fontsize=12, va="center", ha="right")

    ax.text(o_2s_positions[0], y_o_2s - 0.30, r"$a_1$", fontsize=12, ha="center", color="0.35")
    ax.text(o_2px, y_o_2p - 0.32, r"$b_1$", fontsize=12, ha="center", color="0.35")
    ax.text(o_2py, y_o_2p - 0.32, r"$b_2$", fontsize=12, ha="center", color="0.35")
    ax.text(o_2pz, y_o_2p - 0.32, r"$a_1$", fontsize=12, ha="center", color="0.35")

    h_salc_positions = draw_orbital_group(
        ax, x_center=x_h2, y=y_h_salc, number_of_orbitals=2, orbital_spacing=1.8, level_color="tab:purple"
    )
    h_a1_salc = h_salc_positions[0]
    h_b2_salc = h_salc_positions[1]

    ax.text(h_a1_salc, y_h_salc + 0.28, r"$\chi_{a_1}$", fontsize=12, ha="center")
    ax.text(h_b2_salc, y_h_salc + 0.28, r"$\chi_{b_2}$", fontsize=12, ha="center")
    ax.text(h_a1_salc, y_h_salc - 0.35, r"$\frac{1}{\sqrt{2}}(1s_1+1s_2)$", fontsize=12, ha="center", color="0.35")
    ax.text(h_b2_salc, y_h_salc - 0.35, r"$\frac{1}{\sqrt{2}}(1s_1-1s_2)$", fontsize=12, ha="center", color="0.35")

    mo_positions = {}

    for level in H2O_MO_LEVELS:
        positions = draw_orbital_group(
            ax,
            x_center=x_mo,
            y=level["energy"],
            number_of_orbitals=1,
            level_width=0.66,
            level_color=get_level_color(level["type"]),
        )
        mo_positions[level["key"]] = positions

        ax.text(x_mo + 0.55, level["energy"], level["label"], fontsize=13, va="center", ha="left")

    if show_correlation_lines:
        a1_line_color = "0.72"
        b1_line_color = "#d7a76d"
        b2_line_color = "#92b6d5"

        a1_mo_keys = ["2a1", "3a1", "4a1"]

        for mo_key in a1_mo_keys:
            y_mo = get_mo_energy(mo_key)
            x_molecular = mo_positions[mo_key][0]

            connect_from_left(
                ax,
                atomic_x=o_2s_positions[0],
                atomic_y=y_o_2s,
                molecular_x=x_molecular,
                molecular_y=y_mo,
                color=a1_line_color,
            )

            connect_from_left(
                ax, atomic_x=o_2pz, atomic_y=y_o_2p, molecular_x=x_molecular, molecular_y=y_mo, color=a1_line_color
            )

            connect_from_right(
                ax,
                fragment_x=h_a1_salc,
                fragment_y=y_h_salc,
                molecular_x=x_molecular,
                molecular_y=y_mo,
                color=a1_line_color,
            )

        b2_mo_keys = ["1b2", "2b2"]

        for mo_key in b2_mo_keys:
            y_mo = get_mo_energy(mo_key)
            x_molecular = mo_positions[mo_key][0]

            connect_from_left(
                ax, atomic_x=o_2py, atomic_y=y_o_2p, molecular_x=x_molecular, molecular_y=y_mo, color=b2_line_color
            )

            connect_from_right(
                ax,
                fragment_x=h_b2_salc,
                fragment_y=y_h_salc,
                molecular_x=x_molecular,
                molecular_y=y_mo,
                color=b2_line_color,
            )

        connect_from_left(
            ax,
            atomic_x=o_2px,
            atomic_y=y_o_2p,
            molecular_x=mo_positions["1b1"][0],
            molecular_y=get_mo_energy("1b1"),
            color=b1_line_color,
        )

    ax.text(x_o, 8.55, "Oxygen AOs", fontsize=15, ha="center", fontweight="bold", color="tab:green")
    ax.text(x_mo, 8.55, r"H$_2$O MOs", fontsize=15, ha="center", fontweight="bold")
    ax.text(x_h2, 8.55, r"H$_2$ fragment SALCs", fontsize=15, ha="center", fontweight="bold", color="tab:purple")

    ax.annotate(
        "",
        xy=(0, 8.15),
        xytext=(0, 0.45),
        arrowprops=dict(arrowstyle="-|>", linewidth=1.5, color="black"),
    )

    ax.text(-0.38, 4.30, "Energy", fontsize=13, rotation=90, va="center", ha="center")

    ax.axvline(3.7, ymin=0.06, ymax=0.91, color="0.88", linewidth=1.0, linestyle=":")
    ax.axvline(7.5, ymin=0.06, ymax=0.91, color="0.88", linewidth=1.0, linestyle=":")

    ax.set_xlim(-0.5, 12.0)
    ax.set_ylim(0.3, 9.0)

    ax.set_title(r"MO energy diagram of H$_2$O", fontsize=18, pad=18)

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.show()


plot_h2o_mo_diagram()
```
![물 분자의 MO energy diagram](/assets/image-103.png)

표의 맨 오른쪽에 있는 `IE`라는 값은 광전자 분광법(photoelectron spectroscopy, PES)으로 알아낸 실험 값입니다. 분자에 빛을 쪼인 다음에 튀어나오는 전자의 운동 에너지를 측정해서 각 오비탈의 에너지를 계산하는 방법이죠. 광전 효과를 이용한 겁니다.

우리가 계산한 값들과 비교해보면 에너지 값과 순서가 꽤 잘 맞는 것을 볼 수 있습니다. 약간의 오차는 있긴 하지만 $\beta$ 값을 대충 잡고 전자 간의 반발도 무시한 근사적인 계산이라는 걸 생각하면 대단히 성공적인 모델이라고 할 수 있죠.

[[TIP]]
여기서 분자 오비탈의 에너지를 그 오비탈에서 전자를 떼어내는 데 필요한 에너지(이온화 에너지)로 썼습니다. 이 관계를 Koopmans 정리라고 부르는데, 사실 둘이 정확하게 같은 값은 아닙니다. 전자가 빠지면 남은 전자들이 다시 재배치되면서 오비탈 이완이라고 부르는 현상이 일어나기 때문입니다. 이 정리에 대한 더 자세한 내용은 17장에서 얘기해보겠습니다.
[[/TIP]]


## 물 분자의 모양을 설명하기

아까 우리는 물 분자의 결합각을 모른다고 가정하자고 했습니다. (바로 다음에 104.5를 대입하긴 했지만요.) 산소 원자 하나와 수소 원자 두 개에서 출발해서 물 분자를 만들 때, 분자 오비탈을 형성하는 과정에서 다른 각도로 결합하게 되면 에너지가 달라질 겁니다. 결국 가장 안정한 구조를 찾아갈 텐데, 왜 굽은 구조가 가장 안정하게 될까요?  

직접 결합각을 바꿔가면서 오비탈들의 에너지가 어떻게 변하는지 그려봅시다. 가장 에너지가 낮은 구조가 가장 안정한 구조일테니까요.

```python
import matplotlib.pyplot as plt

angles = np.arange(80, 181, 1)
energies = {"2a1": [], "1b2": [], "3a1": [], "1b1": []}
total = []
for theta in angles:
    blocks = water_blocks(theta)
    E_a1 = np.linalg.eigvalsh(blocks["a1"])
    E_b2 = np.linalg.eigvalsh(blocks["b2"])
    E_b1 = np.linalg.eigvalsh(blocks["b1"])
    occ = {"2a1": E_a1[0], "1b2": E_b2[0], "3a1": E_a1[1], "1b1": E_b1[0]}
    for mo in energies:
        energies[mo].append(occ[mo])
    total.append(2 * sum(occ.values()))

fig, axes = plt.subplots(1, 2, figsize=(11, 5))
ax = axes[0]
for mo, col in [("2a1", "magenta"), ("1b2", "seagreen"), ("3a1", "steelblue"), ("1b1", "crimson")]:
    ax.plot(angles, energies[mo], lw=3, color=col, label=mo)
ax.axvline(104.5, color="gray", ls="--", lw=1)
ax.set_xlabel("H-O-H bond angle (degree)")
ax.set_ylabel("orbital energy (eV)")
ax.legend(fontsize=10)
ax.grid(alpha=0.3)

ax = axes[1]
ax.plot(angles, total, lw=3, color="black")
ax.axvline(104.5, color="gray", ls="--", lw=1)
ax.set_xlabel("H-O-H bond angle (degree)")
ax.set_ylabel("total occupied orbital energy (eV)")
ax.grid(alpha=0.3)

plt.tight_layout()
plt.show()
```
![물 분자의 Walsh 도표](/assets/image-104.png)

왼쪽 그래프를 보면 결합각이 180˚에 가까워질수록 $1b_2$ 오비탈의 에너지는 내려가지만 $3a_1$ 오비탈의 에너지가 올라가는 것을 볼 수 있습니다. 위에서 얘기했듯이 사인과 코사인 함수 때문이죠. 나머지 두 오비탈의 에너지는 거의 변화가 없고요.  

오른쪽 그래프는 모든 점유(occupied) 오비탈의 에너지의 합입니다. 95˚ 근처에서 에너지가 가장 낮네요. 단순하게 생각하면 모든 원자가 서로 가장 멀리 떨어진 직선 모양이 제일 안정해보이지만 오비탈 간의 중첩, 특히 $3a_1$ 오비탈의 에너지 안정화 때문에 물 분자는 이렇게 꺾인 형태를 가지게 되는 거죠. 이렇게 결합각에 따라 오비탈 에너지의 변화를 그린 그림을 Walsh 도표라고 부릅니다.  

물론 아까 말한 것처럼 지금 계산은 근사적인 계산이라 정확하게 결합각 104.5˚를 예측하지는 못했지만 분자가 굽어있는 모양이 더 안정하다는 결론은 같습니다.


## 다음 이야기

일반화학에서 물 분자의 구조에 대해서 배울 때 $sp^3$ 혼성 오비탈 네 개가 정사면체 구조를 이루며 그 중 두 개는 수소와 결합하고 나머지 두 개에는 비공유 전자쌍이 있다고 배웠었죠. 이렇게만 들으니 비공유 전자쌍 두 개는 똑같은 것 같습니다.  
그런데 방금 만든 물 분자의 분자 오비탈 그림을 보면 결합에 참여하지 않는 전자가 들어있는 오비탈은 $1b_1$과 $3a_1$이고 둘은 전혀 다른 오비탈입니다! 대칭 종도 다르고 에너지도 다르죠. 이론적으로만 그런 게 아니라 실험적으로도 12.6 eV와 14.7 eV 두 가지가 발견됩니다.  
그러면 일반화학에서 배운 게 틀린 걸까요? 분자 오비탈 이론과 혼성 오비탈 개념은 서로 모순되는 걸까요? 다음 절에서 이 문제를 살펴보겠습니다.


## 확인 문제

1. `water_block`에서 `b_s`와 `b_p`를 바꿔서 계산해보세요. 분자 오비탈의 순서는 어떻게 변하나요? 총 에너지가 최소가 되는 결합각은요?
2. BeH₂는 원자가 전자가 4개뿐이어서 $2a_1$과 $1b_2$ 오비탈까지만 채워집니다. Walsh 도표를 보고 이 분자의 결합각이 얼마일지 예측해보세요. (코드에서 `total`을 구할 때 $3a_1$과 $1b_1$ 오비탈을 빼면 됩니다.)
3. 암모니아 분자의 원자가 오비탈은 어떤 블록들로 나뉠까요? 수소의 $1s$ 오비탈 SALC와 질소의 $2s$, $2p$ 오비탈을 대칭 종별로 묶어보세요.
