# 14.4. 물 분자의 분자 오비탈

(나중에 물 분자 MO diagram 그릴 때 쓸 코드.)
```python
import matplotlib.pyplot as plt

H2O_MO_LEVELS = [
    {
        "key": "2a1",
        "label": r"$2a_1$",
        "energy": 0.6,
        "symmetry": "a1",
        "type": "bonding",
    },
    {
        "key": "1b2",
        "label": r"$1b_2$",
        "energy": 3.2,
        "symmetry": "b2",
        "type": "bonding",
    },
    {
        "key": "3a1",
        "label": r"$3a_1$",
        "energy": 4.5,
        "symmetry": "a1",
        "type": "nonbonding",
    },
    {
        "key": "1b1",
        "label": r"$1b_1$",
        "energy": 4.8,
        "symmetry": "b1",
        "type": "nonbonding",
    },
    {
        "key": "4a1",
        "label": r"$4a_1$",
        "energy": 7.2,
        "symmetry": "a1",
        "type": "antibonding",
    },
    {
        "key": "2b2",
        "label": r"$2b_2$",
        "energy": 7.8,
        "symmetry": "b2",
        "type": "antibonding",
    },
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


def draw_orbital_group(ax, x_center, y, num_orb=1, orb_sp=0.60, level_width=0.46, level_color="black"):

    if num_orb == 1:
        x_positions = [x_center]
    else:
        start = x_center - orb_sp * (num_orb - 1) / 2
        x_positions = [start + i * orb_sp for i in range(num_orb)]

    for x in x_positions:
        draw_energy_line(ax, x_center=x, y=y, width=level_width, color=level_color)

    return x_positions


def draw_correlation_line(ax, x1, y1, x2, y2, color="0.75", linewidth=1.0, linestyle="--"):
    ax.plot([x1, x2], [y1, y2], color=color, linewidth=linewidth, linestyle=linestyle, zorder=1)


def connect_from_left(ax, a_x, a_y, m_x, m_y, color="0.75"):
    draw_correlation_line(ax, x1=a_x + 0.23, y1=a_y, x2=m_x - 0.31, y2=m_y, color=color)


def connect_from_right(ax, frag_x, frag_y, m_x, m_y, color="0.75"):
    draw_correlation_line(ax, x1=frag_x - 0.23, y1=frag_y, x2=m_x + 0.31, y2=m_y, color=color)


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
    x_o = 1.6
    x_mo = 5.0
    x_h2 = 9.0

    y_o_2s = 2.0
    y_o_2p = 4.8
    y_h_s = 4.7
    y_h_s_star = 5.7

    fig, ax = plt.subplots(figsize=(9, 10))

    o_2s_positions = draw_orbital_group(ax, x_center=x_o, y=y_o_2s, num_orb=1, level_color="tab:green")
    o_2p_positions = draw_orbital_group(ax, x_center=x_o, y=y_o_2p, num_orb=3, orb_sp=0.8, level_color="tab:green")

    o_2px = o_2p_positions[0]
    o_2py = o_2p_positions[1]
    o_2pz = o_2p_positions[2]

    ax.text(o_2px, y_o_2p + 0.2, r"$2p_x$", fontsize=12, ha="center")
    ax.text(o_2py, y_o_2p + 0.2, r"$2p_y$", fontsize=12, ha="center")
    ax.text(o_2pz, y_o_2p + 0.2, r"$2p_z$", fontsize=12, ha="center")
    ax.text(x_o, y_o_2s + 0.2, r"$2s$", fontsize=12, ha="center")

    ax.text(o_2s_positions[0], y_o_2s - 0.3, r"$a_1$", fontsize=10, ha="center", color="0.35")
    ax.text(o_2px, y_o_2p - 0.3, r"$b_1$", fontsize=10, ha="center", color="0.35")
    ax.text(o_2py, y_o_2p - 0.3, r"$b_2$", fontsize=10, ha="center", color="0.35")
    ax.text(o_2pz, y_o_2p - 0.3, r"$a_1$", fontsize=10, ha="center", color="0.35")

    h_s_pos = draw_orbital_group(ax, x_center=x_h2, y=y_h_s, num_orb=1, orb_sp=1.5, level_color="tab:purple")
    h_s_star_pos = draw_orbital_group(ax, x_center=x_h2, y=y_h_s_star, num_orb=1, orb_sp=1.5, level_color="tab:purple")

    h_a1_salc = h_s_pos[0]
    h_b2_salc = h_s_star_pos[0]

    ax.text(h_a1_salc, y_h_s + 0.2, r"$\chi_{a_1}$", fontsize=12, ha="center")
    ax.text(h_b2_salc, y_h_s_star + 0.2, r"$\chi_{b_2}$", fontsize=12, ha="center")
    ax.text(h_a1_salc, y_h_s - 0.3, r"$\frac{1}{\sqrt{2}}(1s_1+1s_2)$", fontsize=10, ha="center", color="0.35")
    ax.text(h_b2_salc, y_h_s_star - 0.3, r"$\frac{1}{\sqrt{2}}(1s_1-1s_2)$", fontsize=10, ha="center", color="0.35")

    mo_positions = {}

    for level in H2O_MO_LEVELS:
        positions = draw_orbital_group(
            ax,
            x_center=x_mo,
            y=level["energy"],
            num_orb=1,
            level_width=0.66,
            level_color=get_level_color(level["type"]),
        )

        mo_positions[level["key"]] = positions

        ax.text(x_mo + 0.55, level["energy"], level["label"], fontsize=13, va="center", ha="left")

    if show_correlation_lines:
        a1_mo_keys = ["2a1", "3a1", "4a1"]

        for mo_key in a1_mo_keys:
            y_mo = get_mo_energy(mo_key)
            x_mo = mo_positions[mo_key][0]

            connect_from_left(ax, a_x=o_2s_positions[0], a_y=y_o_2s, m_x=x_mo, m_y=y_mo, color="0.72")
            connect_from_left(ax, a_x=o_2pz, a_y=y_o_2p, m_x=x_mo, m_y=y_mo, color="0.72")
            connect_from_right(ax, frag_x=h_a1_salc, frag_y=y_h_s, m_x=x_mo, m_y=y_mo, color="0.72")

        b2_mo_keys = ["1b2", "2b2"]

        for mo_key in b2_mo_keys:
            y_mo = get_mo_energy(mo_key)
            x_mo = mo_positions[mo_key][0]

            connect_from_left(ax, a_x=o_2py, a_y=y_o_2p, m_x=x_mo, m_y=y_mo, color="0.72")
            connect_from_right(ax, frag_x=h_b2_salc, frag_y=y_h_s_star, m_x=x_mo, m_y=y_mo, color="0.72")

        connect_from_left(ax, a_x=o_2px, a_y=y_o_2p, m_x=mo_positions["1b1"][0], m_y=get_mo_energy("1b1"), color="0.72")

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

    ax.axvline(3.2, ymin=0.06, ymax=0.91, color="0.88", linewidth=1.0, linestyle=":")
    ax.axvline(7.0, ymin=0.06, ymax=0.91, color="0.88", linewidth=1.0, linestyle=":")
    ax.set_xlim(0, 11)
    ax.set_ylim(0, 9.0)
    ax.set_title(r"MO energy diagram of H$_2$O", fontsize=18, pad=18)
    ax.set_xticks([])
    ax.set_yticks([])

    for spine in ax.spines.values():
        spine.set_visible(False)

    plt.tight_layout()
    plt.show()


plot_h2o_mo_diagram(show_correlation_lines=True)
```
