import pandas as pd
import matplotlib.pyplot as plt
import numpy as np

# Read the CSV file (space-separated)
df = pd.read_csv('iv_results.csv', sep='\s+', engine='python')

# Get the time column
time = df.iloc[:, 0]

# Determine the number of I-V pairs
n_pairs = (len(df.columns) - 1) // 2

# Create figure with three subplots
fig, axes = plt.subplots(1, 3, figsize=(15, 5))

# Plot I vs t
for i in range(n_pairs):
    I_col = df.columns[1 + 2*i]  # I columns are at odd indices (1, 3, 5, ...)
    axes[0].plot(time*1000, df[I_col], label=i, marker='o', markersize=3)
axes[0].set_xlabel('Time (ms)')
axes[0].set_ylabel('Current (A)')
axes[0].set_title('Current')
#axes[0].legend()
axes[0].grid(True)

# Plot V vs t
for i in range(n_pairs):
    V_col = df.columns[2 + 2*i]  # V columns are at even indices (2, 4, 6, ...)
    axes[1].plot(time*1000, df[V_col]*1000*2, label=i, marker='o', markersize=3)
axes[1].set_xlabel('Time (ms)')
axes[1].set_ylabel('Voltage (mV)')
axes[1].set_title('Voltage')
#axes[1].legend()
axes[1].grid(True)

# Plot V vs I
for i in range(n_pairs):
    I_col = df.columns[1 + 2*i]
    V_col = df.columns[2 + 2*i]
    axes[2].plot(time*1000, df[I_col]*df[V_col]*2, label=i, marker='o', markersize=3)
axes[2].set_xlabel('Time (ms)')
axes[2].set_ylabel('Losses (W)')
axes[2].set_title('Losses')
#axes[2].legend()
axes[2].grid(True)

plt.tight_layout()
plt.savefig('iv_plots.png', dpi=300, bbox_inches='tight')
plt.show()
