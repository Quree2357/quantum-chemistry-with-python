# 14.3. 이핵 이원자 분자

<a target="_blank" rel="noopener noreferrer" href="https://colab.research.google.com/github/Quree2357/quantum-chemistry-with-python/blob/main/scripts/14-03.ipynb">![Open In Colab](https://colab.research.google.com/assets/colab-badge.svg)</a>

이번에는 서로 다른 두 원자가 결합한 경우를 봅시다. 동핵 이원자 분자에서는 두 원자의 $\alpha$가 같았기 때문에 분자 오비탈이 골고루 섞였습니다. 하지만 이제는 두 원자가 달라서 14.1절에서 본 것처럼 결합성 오비탈과 반결합성 오비탈이 한쪽으로 쏠리게 될 겁니다. 이것이 실제 분자에서 어떤 모습으로 나타나는지 HF와 CO 분자를 통해 살펴보겠습니다.

참고로, 이핵 이원자 분자에서는 반전 중심이 사라져서 분자 오비탈에 아래첨자 $g$와 $u$를 붙일 수 없습니다. 그래서 그냥 번호만 붙이죠.


## 원자는 강약약강이다

14.1절의 확인 문제에서 LiH와 HF 분자를 비교해봤었죠. 수소의 $1s$ 오비탈은 -13.6 eV, 리튬의 $2s$ 오비탈은 -5.39 eV, 그리고 플루오린의 $2p$ 오비탈은 -18.65 eV의 에너지를 가집니다. 결합성 오비탈에서 수소가 차지하는 비중을 계산해보죠. 14.1절의 `orbital_mix` 함수를 그대로 쓰겠습니다. 지금은 $\beta$ 값을 정확히 알 수 없으니 여러 값들을 넣어보겠습니다.

```python
import numpy as np


def orbital_mix(alpha_A, alpha_B, beta):
    H = np.array([[alpha_A, beta], [beta, alpha_B]])
    E, C = np.linalg.eigh(H)
    return E, C


alpha_H, alpha_Li, alpha_F = -13.6, -5.39, -18.65

print(f"{'beta':>5} {'LiH':>7} {'HF':>7}")
for beta in [-3.0, -4.0, -5.0, -6.0, -7.0]:
    _, C_LiH = orbital_mix(alpha_H, alpha_Li, beta)
    _, C_HF = orbital_mix(alpha_H, alpha_F, beta)
    print(f"{beta:5.1f} {C_LiH[0, 0]**2:7.1%} {C_HF[0, 0]**2:7.1%}")
```
```
 beta     LiH      HF
 -3.0   90.4%   17.8%
 -4.0   85.8%   23.3%
 -5.0   81.7%   27.5%
 -6.0   78.2%   30.6%
 -7.0   75.3%   33.0%
```
LiH 분자에서는 결합성 오비탈의 약 75~90%가 수소에 있고, HF 분자에서는 20~30% 정도만 있네요. 결합성 오비탈에 전자 두 개가 들어갈 테니 LiH에서는 수소 쪽에 전자가 쏠리고 HF에서는 반대로 플루오린 쪽에 전자가 쏠릴 겁니다.  
전기음성도를 봐도 같은 경향이 나오죠. 수소의 전기음성도가 2.2인데 플루오린은 4.0이라 수소보다 전자를 더 잘 끌어당기고 리튬의 전기음성도는 0.98이라 수소가 더 전자를 잘 끌어당기니까요. 


## 섞이지 않는 오비탈

HF 분자의 분자 오비탈을 조금 더 자세히 봅시다. 수소의 원자가 오비탈은 $1s$ 하나뿐이고, 플루오린에는 $2s$, $2p_x$, $2p_y$, $2p_z$ 네 개가 있습니다.  

먼저 대칭성부터 살펴보죠. 결합 축을 $z$축으로 두면 14.2절에서 본 것처럼 플루오린의 $2p_x$와 $2p_y$는 결합 축을 포함하는 거울면에 대해 부호가 바뀌는데 수소의 $1s$ 오비탈은 그대로입니다. 그러니 이들은 서로 섞이지 않고 원자 오비탈 그대로 분자 오비탈을 형성합니다. 이런 오비탈을 **비결합 오비탈(nonbonding orbital)**이라고 부르죠.  

다음은 에너지를 살펴볼까요? 플루오린의 $2s$ 오비탈의 에너지는 -40.17 eV로 수소의 $1s$ 오비탈과 26 eV 정도가 차이 납니다. 그래서 이것도 거의 섞이지 않고 사실상 비결합 오비탈로 남죠. 결국 수소의 $1s$ 오비탈과 제대로 섞이는 것은 플루오린의 $2p_z$ 오비탈 하나뿐이고, 이 둘이 결합성 오비탈과 반결합성 오비탈을 만듭니다. 그림으로 그려보면 이렇게 되죠.

(그림 그리고 코드 지우기)
```python
import matplotlib.pyplot as plt

HF_MO_LEVELS = [
    {
        "key": "n_F_2s",
        "label": r"$2\sigma$  (nonbonding)",
        "energy": 1.3,
        "degeneracy": 1,
        "type": "nonbonding",
    },
    {
        "key": "sigma_HF",
        "label": r"$3\sigma$  (bonding)",
        "energy": 3.5,
        "degeneracy": 1,
        "type": "bonding",
    },
    {
        "key": "n_F_2pxy",
        "label": r"$1\pi$  (nonbonding)",
        "energy": 4.6,
        "degeneracy": 2,
        "type": "nonbonding",
    },
    {
        "key": "sigma_star_HF",
        "label": r"$4\sigma^{*}$  (antibonding)",
        "energy": 6.5,
        "degeneracy": 1,
        "type": "antibonding",
    },
]


def distribute_hund(number_of_electrons, number_of_orbitals):

    capacity = 2 * number_of_orbitals
    number_of_electrons = min(number_of_electrons, capacity)

    occupation = [0] * number_of_orbitals

    singly_occupied = min(number_of_electrons, number_of_orbitals)

    for i in range(singly_occupied):
        occupation[i] = 1

    remaining = number_of_electrons - singly_occupied

    for i in range(remaining):
        occupation[i] += 1

    return occupation


def fill_hf_atomic_orbitals():
    return {"F": {"2s": [0], "2p": [0, 0, 0]}, "H": {"1s": [0]}}


def fill_hf_molecular_orbitals():
    total_valence_electrons = 0
    remaining = total_valence_electrons
    occupations = {}

    for level in HF_MO_LEVELS:
        degeneracy = level["degeneracy"]
        capacity = 2 * degeneracy
        electrons_here = min(remaining, capacity)
        occupations[level["key"]] = distribute_hund(electrons_here, degeneracy)
        remaining -= electrons_here

    return occupations


def draw_energy_line(ax, x_center, y, width=0.8, color="black", linewidth=3, zorder=3):
    ax.plot(
        [x_center - width / 2, x_center + width / 2],
        [y, y],
        color=color,
        linewidth=linewidth,
        solid_capstyle="butt",
        zorder=zorder,
    )


def draw_electron_arrow(ax, x, y, spin="up", color="tab:red", height=0.42):
    if spin == "up":
        y_start = y + 0.03
        y_end = y + height
    else:
        y_start = y + height
        y_end = y + 0.03

    ax.annotate(
        "",
        xy=(x, y_end),
        xytext=(x, y_start),
        arrowprops=dict(
            arrowstyle="-|>",
            color=color,
            linewidth=1.35,
            mutation_scale=10,
        ),
        zorder=5,
    )


def draw_orbital_group(ax, x_center, y, occupations, orbital_spacing=0.58, level_width=0.46, level_color="black"):
    number_of_orbitals = len(occupations)

    if number_of_orbitals == 1:
        x_positions = [x_center]
    else:
        start = x_center - orbital_spacing * (number_of_orbitals - 1) / 2

        x_positions = [start + i * orbital_spacing for i in range(number_of_orbitals)]

    for x in x_positions:
        draw_energy_line(ax, x_center=x, y=y, width=level_width, color=level_color)

    return x_positions


def draw_correlation_line(ax, x1, y1, x2, y2, color="0.72", linewidth=1.1, linestyle="--"):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle=linestyle, zorder=1)


def find_mo_energy(key):
    for level in HF_MO_LEVELS:
        if level["key"] == key:
            return level["energy"]

    raise KeyError(f"Unknown MO key: {key}")


def level_color(level_type):
    colors = {
        "bonding": "tab:blue",
        "nonbonding": "black",
        "antibonding": "tab:red",
    }

    return colors[level_type]


def plot_hf_mo_diagram(show_correlation_lines=True):
    atomic_occupation = fill_hf_atomic_orbitals()
    mo_occupation = fill_hf_molecular_orbitals()

    x_f = 1.5
    x_mo = 5.0
    x_h = 8.5

    # 원자 오비탈 에너지
    y_f_2s = 1.3
    y_f_2p = 4.6
    y_h_1s = 5.3

    fig, ax = plt.subplots(figsize=(8, 10))

    f_2s_x = draw_orbital_group(
        ax,
        x_center=x_f,
        y=y_f_2s,
        occupations=atomic_occupation["F"]["2s"],
        level_color="tab:green",
    )

    f_2p_x = draw_orbital_group(
        ax,
        x_center=x_f,
        y=y_f_2p,
        occupations=atomic_occupation["F"]["2p"],
        orbital_spacing=0.58,
        level_color="tab:green",
    )

    h_1s_x = draw_orbital_group(
        ax,
        x_center=x_h,
        y=y_h_1s,
        occupations=atomic_occupation["H"]["1s"],
        level_color="tab:green",
    )

    mo_positions = {}

    for level in HF_MO_LEVELS:
        positions = draw_orbital_group(
            ax,
            x_center=x_mo,
            y=level["energy"],
            occupations=mo_occupation[level["key"]],
            orbital_spacing=0.72,
            level_width=0.62,
            level_color=level_color(level["type"]),
        )

        mo_positions[level["key"]] = positions

        ax.text(x_mo + 0.95, level["energy"], level["label"], fontsize=12, va="center", ha="left")

    if show_correlation_lines:
        y_mo = find_mo_energy("n_F_2s")

        draw_correlation_line(ax, f_2s_x[0] + 0.23, y_f_2s, mo_positions["n_F_2s"][0] - 0.31, y_mo)

        f_2pz_position = f_2p_x[2]

        for mo_key in ["sigma_HF", "sigma_star_HF"]:
            y_mo = find_mo_energy(mo_key)
            mo_x = mo_positions[mo_key][0]

            draw_correlation_line(ax, f_2pz_position + 0.23, y_f_2p, mo_x - 0.31, y_mo)

            draw_correlation_line(ax, h_1s_x[0] - 0.23, y_h_1s, mo_x + 0.31, y_mo)

        y_pi = find_mo_energy("n_F_2pxy")

        for atomic_x, molecular_x in zip(f_2p_x[:2], mo_positions["n_F_2pxy"]):
            draw_correlation_line(ax, atomic_x + 0.23, y_f_2p, molecular_x - 0.31, y_pi)

    ax.text(x_f, 8.15, "Fluorine", fontsize=15, ha="center", fontweight="bold")
    ax.text(x_mo, 8.15, "HF molecular orbitals", fontsize=15, ha="center", fontweight="bold")
    ax.text(x_h, 8.15, "Hydrogen", fontsize=15, ha="center", fontweight="bold")

    ax.text(x_f - 0.65, y_f_2s, r"$2s$", fontsize=13, va="center", ha="right")
    ax.text(x_f - 1, y_f_2p, r"$2p$", fontsize=13, va="center", ha="right")
    ax.text(x_h + 0.65, y_h_1s, r"$1s$", fontsize=13, va="center", ha="left")

    ax.annotate(
        "",
        xy=(0, 7.8),
        xytext=(0, 0.5),
        arrowprops=dict(arrowstyle="-|>", linewidth=1.5, color="black"),
    )

    ax.text(-0.35, 4.15, "Energy", fontsize=13, rotation=90, va="center", ha="center")

    # 영역 구분선
    ax.axvline(3.0, 0.06, 0.89, color="0.88", linewidth=1.0, linestyle=":")
    ax.axvline(7.0, 0.06, 0.89, color="0.88", linewidth=1.0, linestyle=":")

    ax.set_xlim(-0.5, 10)
    ax.set_ylim(0.5, 8.6)
    ax.set_title("MO energy diagram of HF", fontsize=18, pad=16)
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.show()


plot_hf_mo_diagram()
```
![HF의 MO energy diagram](/assets/image-101.png)

원자가 전자 8개를 아래에서부터 채워나가면 전자 배치는 $(2\sigma)^2(3\sigma)^2(1\pi)^4$가 됩니다. 결합성 오비탈에 전자 두 개가 있으니 결합 차수는 1이죠. 그리고 나머지 전자 여섯 개는 비결합 오비탈에 들어가 있습니다.

Lewis 구조로 HF 분자를 그리면 플루오린 원자에 비공유 전자쌍이 세 개 있죠. 분자 오비탈에서도 비결합 오비탈에 전자 여섯 개가 있습니다. 하지만 세 쌍이 전부 동일한 것은 아닙니다. 그림에서 볼 수 있듯이 $2\sigma$에 들어 있는 한 쌍은 에너지가 낮고, $1\pi$에 들어 있는 나머지 두 쌍은 에너지가 높죠. 실제로 HF 분자가 이온화될 때 떨어지는 전자는 $1\pi$의 전자입니다.  


## 청개구리 같은 CO

CO(일산화탄소) 분자는 원자가 전자가 10개로 질소 분자와 같습니다. 이런 관계를 등전자(isoelectronic)라고 부르죠. 그래서 분자 오비탈도 질소 분자와 비슷하지만 산소의 전기음성도가 더 커서 산소의 원자 오비탈이 더 아래로 내려가 있습니다.

(그림 그리고 코드 지우기)
```python
import matplotlib.pyplot as plt

CO_MO_LEVELS = [
    {
        "key": "3sigma",
        "label": r"$3\sigma$",
        "energy": 1.1,
        "degeneracy": 1,
        "type": "bonding",
    },
    {
        "key": "4sigma",
        "label": r"$4\sigma$",
        "energy": 3.3,
        "degeneracy": 1,
        "type": "antibonding",
    },
    {
        "key": "1pi",
        "label": r"$1\pi$",
        "energy": 4.1,
        "degeneracy": 2,
        "type": "bonding",
    },
    {
        "key": "5sigma",
        "label": r"$5\sigma$",
        "energy": 5.1,
        "degeneracy": 1,
        "type": "bonding",
    },
    {
        "key": "2pi_star",
        "label": r"$2\pi^{*}$",
        "energy": 7.5,
        "degeneracy": 2,
        "type": "antibonding",
    },
    {
        "key": "6sigma_star",
        "label": r"$6\sigma^{*}$",
        "energy": 8.0,
        "degeneracy": 1,
        "type": "antibonding",
    },
]

CO_MO_OCCUPATION = {"3sigma": 2, "4sigma": 2, "1pi": 4, "5sigma": 2, "2pi_star": 0, "6sigma_star": 0}


def draw_energy_line(ax, x_center, y, width=0.8, color="black", linewidth=3, zorder=3):
    ax.plot(
        [x_center - width / 2, x_center + width / 2],
        [y, y],
        color=color,
        linewidth=linewidth,
        solid_capstyle="butt",
        zorder=zorder,
    )


def draw_orbital_group(ax, x_center, y, degeneracy=1, orbital_spacing=0.60, level_width=0.46, level_color="black"):
    if degeneracy == 1:
        x_positions = [x_center]
    else:
        start = x_center - orbital_spacing * (degeneracy - 1) / 2
        x_positions = [start + i * orbital_spacing for i in range(degeneracy)]

    for x in x_positions:
        draw_energy_line(ax, x_center=x, y=y, width=level_width, color=level_color)

    return x_positions


def draw_correlation_line(ax, x1, y1, x2, y2, color="0.72", linewidth=1.0, linestyle="--"):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle=linestyle, zorder=1)


def get_mo_level(key):
    for level in CO_MO_LEVELS:
        if level["key"] == key:
            return level

    raise KeyError(f"Unknown MO key: {key}")


def get_mo_energy(key):
    return get_mo_level(key)["energy"]


def get_level_color(level_type):
    colors = {"bonding": "tab:blue", "nonbonding": "black", "antibonding": "tab:red"}

    return colors[level_type]


def connect_single_ao_to_mo(ax, atomic_x, atomic_y, molecular_x, molecular_y, side="left"):
    if side == "left":
        x1 = atomic_x + 0.23
        x2 = molecular_x - 0.31
    else:
        x1 = atomic_x - 0.23
        x2 = molecular_x + 0.31

    draw_correlation_line(ax, x1=x1, y1=atomic_y, x2=x2, y2=molecular_y)


def plot_co_mo_diagram(show_correlation_lines=True):

    x_c = 1.5
    x_mo = 5.5
    x_o = 9.5

    y_c_2s = 3.0
    y_c_2p = 7.0

    y_o_2s = 1.3
    y_o_2p = 4.5

    fig, ax = plt.subplots(figsize=(8, 10))

    c_2s_x = draw_orbital_group(ax, x_center=x_c, y=y_c_2s, degeneracy=1, level_color="tab:green")
    c_2p_x = draw_orbital_group(ax, x_center=x_c, y=y_c_2p, degeneracy=3, orbital_spacing=0.58, level_color="tab:green")
    c_2px_position = c_2p_x[0]
    c_2py_position = c_2p_x[1]
    c_2pz_position = c_2p_x[2]

    o_2s_x = draw_orbital_group(ax, x_center=x_o, y=y_o_2s, degeneracy=1, level_color="tab:green")
    o_2p_x = draw_orbital_group(ax, x_center=x_o, y=y_o_2p, degeneracy=3, orbital_spacing=0.58, level_color="tab:green")
    o_2px_position = o_2p_x[0]
    o_2py_position = o_2p_x[1]
    o_2pz_position = o_2p_x[2]

    mo_positions = {}

    for level in CO_MO_LEVELS:
        positions = draw_orbital_group(
            ax,
            x_center=x_mo,
            y=level["energy"],
            degeneracy=level["degeneracy"],
            orbital_spacing=0.72,
            level_width=0.62,
            level_color=get_level_color(level["type"]),
        )

        mo_positions[level["key"]] = positions

        ax.text(x_mo + 1, level["energy"], level["label"], fontsize=12, va="center", ha="left")

    if show_correlation_lines:
        for mo_key in ["3sigma", "4sigma"]:
            y_mo = get_mo_energy(mo_key)
            x_mo = mo_positions[mo_key][0]

            connect_single_ao_to_mo(
                ax, atomic_x=c_2s_x[0], atomic_y=y_c_2s, molecular_x=x_mo, molecular_y=y_mo, side="left"
            )
            connect_single_ao_to_mo(
                ax, atomic_x=o_2s_x[0], atomic_y=y_o_2s, molecular_x=x_mo, molecular_y=y_mo, side="right"
            )

        c_pi_ao = c_2pz_position
        o_pi_ao = o_2px_position
        y_5s = get_mo_energy("5sigma")
        x_5s = mo_positions["5sigma"][0]
        y_1p = get_mo_energy("1pi")
        x_1p = mo_positions["1pi"][0]
        x_1p_2 = mo_positions["1pi"][1]
        y_2p = get_mo_energy("2pi_star")
        x_2p = mo_positions["2pi_star"][0]
        x_2p_2 = mo_positions["2pi_star"][1]
        y_6s = get_mo_energy("6sigma_star")
        x_6s = mo_positions["6sigma_star"][0]
        connect_single_ao_to_mo(ax, atomic_x=c_pi_ao, atomic_y=y_c_2p, molecular_x=x_1p, molecular_y=y_1p, side="left")
        connect_single_ao_to_mo(ax, atomic_x=o_pi_ao, atomic_y=y_o_2p, molecular_x=x_1p_2, molecular_y=y_1p, side="right")
        connect_single_ao_to_mo(ax, atomic_x=c_pi_ao, atomic_y=y_c_2p, molecular_x=x_2p, molecular_y=y_2p, side="left")
        connect_single_ao_to_mo(ax, atomic_x=o_pi_ao, atomic_y=y_o_2p, molecular_x=x_2p_2, molecular_y=y_2p, side="right")
        connect_single_ao_to_mo(ax, atomic_x=c_pi_ao, atomic_y=y_c_2p, molecular_x=x_5s, molecular_y=y_5s, side="left")
        connect_single_ao_to_mo(ax, atomic_x=o_pi_ao, atomic_y=y_o_2p, molecular_x=x_5s, molecular_y=y_5s, side="right")
        connect_single_ao_to_mo(ax, atomic_x=c_pi_ao, atomic_y=y_c_2p, molecular_x=x_6s, molecular_y=y_6s, side="left")
        connect_single_ao_to_mo(ax, atomic_x=o_pi_ao, atomic_y=y_o_2p, molecular_x=x_6s, molecular_y=y_6s, side="right")

    ax.text(x_c, 8.7, "Carbon", fontsize=15, ha="center", fontweight="bold", color="black")
    ax.text(x_mo, 8.7, "CO molecular orbitals", fontsize=15, ha="center", fontweight="bold")
    ax.text(x_o, 8.7, "Oxygen", fontsize=15, ha="center", fontweight="bold", color="black")

    ax.text(x_c - 0.65, y_c_2s, r"$2s$", fontsize=13, va="center", ha="right")
    ax.text(x_c - 1.05, y_c_2p, r"$2p$", fontsize=13, va="center", ha="right")

    ax.text(x_o + 0.65, y_o_2s, r"$2s$", fontsize=13, va="center", ha="left")
    ax.text(x_o + 1.05, y_o_2p, r"$2p$", fontsize=13, va="center", ha="left")

    ax.annotate(
        "",
        xy=(0, 8.5),
        xytext=(0, 0.6),
        arrowprops=dict(
            arrowstyle="-|>",
            linewidth=1.5,
            color="black",
        ),
    )

    ax.text(-0.38, 4.3, "Energy", fontsize=13, rotation=90, va="center", ha="center")

    # 영역 구분선
    ax.axvline(3.5, ymin=0.05, ymax=0.90, color="0.88", linewidth=1.0, linestyle=":")
    ax.axvline(8.0, ymin=0.05, ymax=0.90, color="0.88", linewidth=1.0, linestyle=":")

    # 전체 그림 설정
    ax.set_xlim(-0.5, 11.2)
    ax.set_ylim(0.5, 9.0)

    ax.set_title("MO energy diagram of CO", fontsize=18, pad=18)

    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.show()

plot_co_mo_diagram()
```
![CO의 MO energy diagram](/assets/image-102.png)

마찬가지로 각 오비탈에서 탄소와 산소의 비중이 얼마인지 계산해보겠습니다. 이번에는 조금 더 복잡해지는데요. $\beta$값이 음수가 되도록 두 원자의 $2p_z$ 오비탈을 서로를 향하게 잡겠습니다. $\sigma$ 블록에는 탄소의 $2s$와 $2p_z$ 오비탈, 그리고 산소의 $2s$와 $2p_z$ 오비탈 네 개가 들어가서 4×4 행렬이 되고 $\pi$ 블록은 나머지 오비탈들로 이루어진 x, y 각각 하나씩 2×2 행렬 2개가 됩니다. $\beta$ 값은 이번에도 대략적으로 넣어보죠.
```python
a = {"C2s": -19.43, "C2p": -10.66, "O2s": -32.38, "O2p": -15.85}
b_ss, b_sp, b_pp, b_pi = -3, -4, -5, -2.5

H_sigma = np.diag([a["C2s"], a["C2p"], a["O2s"], a["O2p"]])
H_sigma[0, 2], H_sigma[2, 0] = b_ss, b_ss
H_sigma[0, 3], H_sigma[3, 0] = b_sp, b_sp
H_sigma[1, 2], H_sigma[2, 1] = b_sp, b_sp
H_sigma[1, 3], H_sigma[3, 1] = b_pp, b_pp
E_s, C_s = np.linalg.eigh(H_sigma)

E_p, C_p = orbital_mix(a["C2p"], a["O2p"], b_pi)

levels = []
for i in range(4):
    levels.append((E_s[i], f"{i+3}σ", 2, C_s[0, i] ** 2 + C_s[1, i] ** 2))
for i in range(2):
    levels.append((E_p[i], f"{i+1}π", 4, C_p[0, i] ** 2))
levels.sort()

n = 10
print(f"{'MO':>4} {'E (eV)':>8} {'C weight':>10} {'e':>3}")
for E, name, cap, w_C in levels:
    k = min(n, cap)
    n -= k
    print(f"{name:>4} {E:8.2f} {w_C:10.1%} {k:3d}")
```
```
  MO   E (eV)   C weight   e
  3σ   -33.87       8.5%   2
  4σ   -21.83      53.4%   2
  1π   -16.86      14.0%   4
  5σ   -15.96      68.3%   2
  2π    -9.65      86.0%   0
  6σ    -6.66      69.8%   0
```
`C weight`는 각 분자 오비탈에서 탄소가 차지하는 비중입니다. 결합성 $\pi$ 오비탈인 $1\pi$에는 산소가 86%를 기여하고 있습니다. 산소의 오비탈 에너지가 더 낮으니 그럴듯해 보이죠. 그런데 정작 HOMO인 $5\sigma$ 오비탈은 약 68%의 비중을 탄소가 차지합니다! 전기음성도가 더 큰 산소 쪽으로 전자가 쏠려 있을 것 같은데 오히려 탄소 쪽으로 쏠려 있는 거죠. ($3\sigma$와 $4\sigma$의 경우는 이 모델이 잘 안 맞아서 그냥 넘어가겠습니다.) 

$5\sigma$ 분자 오비탈을 좀 더 자세히 봅시다.
```python
c = C_s[:, 2] * np.sign(C_s[0, 2])
print(f"5σ: C 2s {c[0]:+.2f}, C 2pz {c[1]:+.2f}, O 2s {c[2]:+.2f}, O 2p {c[3]:+.2f}")
```
```
5σ: C 2s +0.64, C 2pz -0.52, O 2s +0.01, O 2p -0.56
```
탄소의 $2s$의 부호는 플러스인데 탄소의 $2p_z$와 산소의 $2p_z$ 오비탈의 부호는 마이너스죠. 이렇게 섞이면 탄소 원자의 $2p_z$ 오비탈의 두 로브 중에서 산소 원자의 반대쪽에 있는 것이 더 커지게 됩니다. 즉, 탄소 원자에 있는 비공유 전자쌍에 해당하는 전자가 여기에 들어있는 거죠. 이렇게 $s$와 $p$가 섞여서 만들어지는 오비탈, 어디서 많이 들어본 얘기 아닌가요? 14.5절에서 더 자세히 얘기해보도록 하겠습니다.

또 이 전자쌍은 CO 분자의 또 다른 성질을 만들어냅니다. CO는 금속과 잘 결합하는데, 특히 혈액 속 헤모글로빈의 철 원자에 결합해서 산소 운반을 막는 물질로 잘 알려져 있죠. CO의 HOMO 오비탈에 있는 전자가 탄소 쪽으로 쏠려있기 때문에 CO 분자가 금속 원자에 결합할 때는 전기음성도가 더 큰 산소가 아닌 탄소 쪽으로 결합하게 됩니다. 이렇게 해야 금속의 오비탈과 더 많이 겹쳐질 수 있거든요.

[[TIP]]
CO 분자는 뒤집어진 쌍극자 모멘트로도 유명합니다. 산소의 전기음성도가 더 크니 $\text{C}^{\delta +}\text{O}^{\delta -}$가 되어야 할 것 같지만 실제로는 반대로 $\text{C}^{\delta -}\text{O}^{\delta +}$이고 그 크기도 0.11 D 정도로 아주 작습니다. 단순한 계산으로는 이 현상을 전혀 예측하지 못하죠. 전자들이 서로를 피해다니는 전자 상관 효과까지 고려해서 계산해야 실험과 같은 결과가 나옵니다. 전기음성도만으로 전자의 거동을 판단하면 안 된다는 좋은 예시입니다.
[[/TIP]]


## 다음 이야기

이원자 분자의 경우를 동핵일 때와 이핵일 때 모두 살펴봤습니다. 이제는 원자가 3개인 분자로 넘어가봅시다. 물 분자는 어떨까요? 13.4절에서 물 분자의 행렬이 세 개의 블록으로 쪼개진다는 것을 봤었죠. 각 블록이 어떤 분자 오비탈이 되는지 확인해봅시다.


## 확인 문제

1. N₂ 분자의 $1\pi_u$는 두 원자에 똑같이 퍼져 있는데 CO 분자의 $1\pi$는 86%가 산소에서 기여합니다. 같은 등전자 분자인데 왜 이런 차이가 생길까요?
2. 시안화 이온 $\text{CN}^-$도 CO와 등전자 분자입니다. 시안화 이온이 금속에 결합할 때 어느 쪽으로 결합할지 예측해보세요.
3. CO에서 전자 하나를 떼어내 $\text{CO}^+$를 만들면 결합 길이가 짧아집니다. 위 계산 결과에서 $5\sigma$의 계수 부호를 보고 C의 $2s$와 O의 $2p_z$, C의 $2p_z$와 O의 $2p_z$ 쌍이 각각 결합성인지 반결합성인지 따져보세요. 이 결과와 연관지어 $\text{CO}^+$의 결합 길이를 설명해보세요.
