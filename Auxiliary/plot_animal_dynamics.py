"""
Plot AFT population dynamics from ORNL_animal_test_Animal_yearly.csv.
Also plots plant C from ORNL_animal_test_Ecosystem_yearly.csv.

Panel layout (2 × 2):
  [0,0]  Population density        — Deer (left axis) + Wolf (right axis)
  [0,1]  Ecosystem plant C         — from Ecosystem_yearly.csv
  [1,0]  Annual intake             — Deer plant intake + Wolf prey intake combined
  [1,1]  Starvation mortality rate — both AFTs
"""

import pandas as pd
import matplotlib
matplotlib.use('Agg')   # non-interactive; remove to pop up a window
import matplotlib.pyplot as plt
import sys
import os

# ── data paths ────────────────────────────────────────────────────────────────
script_dir  = os.path.dirname(os.path.abspath(__file__))
output_dir  = os.path.join(script_dir, '..', 'output')
run_id      = 'ORNL_animal_test'

default_ani = os.path.join(output_dir, f'{run_id}_Animal_yearly.csv')
default_eco = os.path.join(output_dir, f'{run_id}_Ecosystem_yearly.csv')

ani_file = sys.argv[1] if len(sys.argv) > 1 else default_ani
eco_file = sys.argv[2] if len(sys.argv) > 2 else default_eco

# ── load animal data ──────────────────────────────────────────────────────────
df = pd.read_csv(ani_file, skipinitialspace=True)
df.columns = df.columns.str.strip()

deer = df[df['AFT'] == 0].sort_values('year')
wolf = df[df['AFT'] == 1].sort_values('year')

# ── load ecosystem data ───────────────────────────────────────────────────────
eco = pd.read_csv(eco_file, skipinitialspace=True)
eco.columns = eco.columns.str.strip()
eco = eco.sort_values('year')

# ── style ─────────────────────────────────────────────────────────────────────
COLORS = {0: '#2196F3', 1: '#F44336'}   # blue = deer, red = wolf
COLOR_ECO = '#4CAF50'                   # green for ecosystem

fig, axes = plt.subplots(2, 2, figsize=(12, 8), sharex=True)
fig.suptitle('Animal population dynamics — ORNL animal test', fontsize=13)

# ── [0,0]  Population density: dual y-axis ────────────────────────────────────
ax_deer = axes[0, 0]
ax_wolf = ax_deer.twinx()

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

# ── [0,1]  Ecosystem plant C ──────────────────────────────────────────────────
ax_pc = axes[0, 1]
ax_pc.plot(eco['year'], eco['plantC'], color=COLOR_ECO, linewidth=1.8)
ax_pc.set_ylabel('Plant C\n(kg C m⁻²)', fontsize=9)
ax_pc.set_title('Plant carbon', fontsize=9)
ax_pc.grid(True, linestyle='--', alpha=0.4)

# ── [1,0]  Combined intake: Deer plant intake + Wolf prey intake ──────────────
ax_int = axes[1, 0]

ln_d, = ax_int.plot(deer['year'], deer['IntakePlant'],
                    color=COLORS[0], linewidth=1.8, label='Deer — plant intake')
ln_w, = ax_int.plot(wolf['year'], wolf['IntakePrey'],
                    color=COLORS[1], linewidth=1.8, linestyle='--',
                    label='Wolf — prey intake')

ax_int.set_ylabel('Annual intake\n(kg C ind⁻¹ yr⁻¹)', fontsize=9)
ax_int.set_title('Annual intake (deer: plant; wolf: prey)', fontsize=9)
ax_int.grid(True, linestyle='--', alpha=0.4)
ax_int.legend(fontsize=8)

# ── [1,1]  Starvation mortality ───────────────────────────────────────────────
ax_mu = axes[1, 1]
for aft_id, sub, label in [(0, deer, 'Deer'), (1, wolf, 'Wolf')]:
    ax_mu.plot(sub['year'], sub['mu_starve'],
               color=COLORS[aft_id], linewidth=1.8,
               linestyle='-' if aft_id == 0 else '--', label=label)
ax_mu.set_ylabel('Starvation mortality\n(day⁻¹)', fontsize=9)
ax_mu.set_title('Starvation mortality rate', fontsize=9)
ax_mu.grid(True, linestyle='--', alpha=0.4)
ax_mu.legend(fontsize=8)

# ── x labels on bottom row only ───────────────────────────────────────────────
for ax in axes[1]:
    ax.set_xlabel('Year', fontsize=9)

# ── save ──────────────────────────────────────────────────────────────────────
plt.tight_layout()
out_png = ani_file.replace('_Animal_yearly.csv', '_dynamics.png')
plt.savefig(out_png, dpi=150)
print(f'Saved: {out_png}')
plt.show()
