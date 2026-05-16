"""
Plot AFT population dynamics from ORNL_animal_test_Animal_yearly.csv.
Panels: population density, plant intake, prey intake, starvation mortality.
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')   # non-interactive; remove to pop up a window
import matplotlib.pyplot as plt
import sys
import os

# ── data path ────────────────────────────────────────────────────────────────
script_dir = os.path.dirname(os.path.abspath(__file__))
default_csv = os.path.join(script_dir, '..', 'output',
                           'ORNL_animal_test_Animal_yearly.csv')
csv_file = sys.argv[1] if len(sys.argv) > 1 else default_csv

# ── load ─────────────────────────────────────────────────────────────────────
df = pd.read_csv(csv_file, skipinitialspace=True)
df.columns = df.columns.str.strip()

AFT_LABELS = {0: 'Deer (AFT 0, herbivore)', 1: 'Wolf (AFT 1, carnivore)'}
COLORS     = {0: '#2196F3', 1: '#F44336'}   # blue, red

# ── figure ───────────────────────────────────────────────────────────────────
fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
fig.suptitle('Animal population dynamics — ORNL animal test', fontsize=13)

panel_cfg = [
    # (ax,        column,         ylabel,                       scale, unit)
    (axes[0, 1], 'IntakePlant', 'Annual plant intake',         1e0,  'kg DM m⁻² yr⁻¹'),
    (axes[1, 0], 'IntakePrey',  'Annual prey intake',          1e0,  'kg C m⁻² yr⁻¹'),
    (axes[1, 1], 'mu_starve',   'Starvation mortality rate',   1e0,  'day⁻¹'),
]

# ── population density panel: dual y-axis ────────────────────────────────────
ax_deer = axes[0, 0]
ax_wolf = ax_deer.twinx()

deer = df[df['AFT'] == 0].sort_values('year')
wolf = df[df['AFT'] == 1].sort_values('year')

ln1, = ax_deer.plot(deer['year'], deer['nindivs'] * 1e4,
                    color=COLORS[0], linewidth=1.8, label='Deer (left)')
ln2, = ax_wolf.plot(wolf['year'], wolf['nindivs'] * 1e4,
                    color=COLORS[1], linewidth=1.8, linestyle='--', label='Wolf (right)')

ax_deer.set_ylabel('Deer density  (ind ha⁻¹)', color=COLORS[0], fontsize=9)
ax_wolf.set_ylabel('Wolf density  (ind ha⁻¹)', color=COLORS[1], fontsize=9)
ax_deer.tick_params(axis='y', labelcolor=COLORS[0])
ax_wolf.tick_params(axis='y', labelcolor=COLORS[1])
ax_deer.set_title('Population density', fontsize=9)
ax_deer.grid(True, linestyle='--', alpha=0.4)
ax_deer.legend(handles=[ln1, ln2], fontsize=8, loc='upper right')

# ── remaining panels ─────────────────────────────────────────────────────────
for ax, col, ylabel, scale, unit in panel_cfg:
    for aft_id, label in AFT_LABELS.items():
        sub = df[df['AFT'] == aft_id].sort_values('year')
        if sub.empty:
            continue
        ax.plot(sub['year'], sub[col] * scale,
                color=COLORS[aft_id], label=label, linewidth=1.8)
    ax.set_ylabel(f'{ylabel}\n({unit})', fontsize=9)
    ax.grid(True, linestyle='--', alpha=0.4)
    ax.legend(fontsize=8)

for ax in axes[1]:
    ax.set_xlabel('Year', fontsize=9)

plt.tight_layout()
out_png = csv_file.replace('.csv', '_dynamics.png')
plt.savefig(out_png, dpi=150)
print(f'Saved: {out_png}')
plt.show()
