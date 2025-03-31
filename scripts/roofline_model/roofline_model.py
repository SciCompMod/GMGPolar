import matplotlib.pyplot as plt
import numpy as np
from matplotlib.lines import Line2D

### Data ###
title = "Roofline Model - Czarny Geometry"
x_axis_label = "Operational Intensity [FLOPS/Byte]"
y_axis_label = "Performance [GFLOPS]"
# Measured by Likwid
# maximum_performance = 1256.6 * 1000  # Maximum performance in MFLOP/s
# maximum_bandwidth = 355.6 * 1000  # Maximum bandwidth in MByte/s
# Given in cluster documentation
maximum_performance = 2.6 * 32 * 56 * 1000
maximum_bandwidth = 281.0 * 1000

### Start Plotting ###
x_min = 1/8
x_max = 128
x_ticks = [1/8, 1/4, 1/2, 1, 2, 4, 8, 16, 32, 64, 128]
x_tick_labels = ["1/8", "1/4", "1/2", "1", "2", "4", "8", "16", "32", "64", "128"]

y_min = 1
y_max = 16384
y_ticks = [1, 4, 16, 64, 256, 1024, 4096, 16384]
y_tick_labels = [f"{int(tick)}" for tick in y_ticks]

# Convert maximum performance to GFLOPS and maximum bandwidth to GB/s
max_perf_gflops = maximum_performance / 1000  # GFLOPS
max_bw_gbs = maximum_bandwidth / 1000  # GB/s

# Operational intensity (OI) range for the plot
oi = np.logspace(np.log10(x_min), np.log10(x_max), 500)

# Roofline calculations
# Memory-bound (Bandwidth roofline)
memory_bound = max_bw_gbs * oi

# Performance-bound (Compute roofline)
performance_bound = np.full_like(oi, max_perf_gflops)

# Intersection point between memory and compute roofs
intersection_oi = max_perf_gflops / max_bw_gbs

# Plot
plt.figure(figsize=(6, 7))

# Memory Bound
x1 = np.array([0, intersection_oi])
y1 = np.array([0, max_perf_gflops])
plt.loglog(x1, y1, 'k-', label="_nolegend_")
# Computation Bound
x1 = np.array([intersection_oi, x_max])
y1 = np.array([max_perf_gflops, max_perf_gflops])
plt.loglog(x1, y1, 'k-', label="_nolegend_")

### No SIMD line ###
no_simd = 2.6 * 56
intersection_no_simd = no_simd / max_bw_gbs
x1 = np.array([intersection_no_simd, x_max])
y1 = np.array([no_simd, no_simd])
plt.loglog(x1, y1, 'g-', label="_nolegend_")  # Green line

# Add text above the NO SIMD line
plt.text(
    intersection_no_simd + 20,  # X position, slightly to the right of the intersection
    no_simd * 1.1,                # Y position, a bit above the line
    'no SIMD',                    # Text label
    color='green',                 # Text color (match the line color)
    fontsize=12,                   # Font size
    verticalalignment='bottom',    # Vertical alignment (place text above the line)
    horizontalalignment='left'    # Horizontal alignment (adjust to the left)
)


# Serial Line in Purple
serial = 2.6
intersection_serial = serial / max_bw_gbs
x1 = np.array([intersection_serial, x_max])
y1 = np.array([serial, serial])
plt.loglog(x1, y1, color='purple', linestyle='-', label="_nolegend_")

# Add text above the NO SIMD line (Purple)
plt.text(
    intersection_serial + 25,  # X position, slightly to the right of the intersection
    serial * 1.1,                # Y position, a bit above the line
    'serial',                    # Text label
    color='purple',                # Text color (match the line color)
    fontsize=12,                   # Font size
    verticalalignment='bottom',    # Vertical alignment (place text above the line)
    horizontalalignment='left'    # Horizontal alignment (adjust to the left)
)


data = [
    # [op_ins, app_perf, color, shape]
    [57.6606, 32448.0051, 'g', 'o'],  # Give (769 x 1024)
    [26.4743, 51214.0858, 'b', 'o'],  # Give (1537 x 2048)
    [9.4655, 67896.1025, 'orange', 'o'],  # Give (3073 x 4096)
    [6.3435, 79118.8143, 'r', 'o'],  # Give (6145 x 8192)
    [32.5134, 16787.3805, 'g', '^'],  # Take (769 x 1024)
    [4.9066, 28809.3174, 'b', '^'],  # Take (1537 x 2048)
    [3.2515, 45878.1637, 'orange', '^'],  # Take (3073 x 4096)
    [2.9040, 54252.0375, 'r', '^'],  # Take (6145 x 8192)
]

for op_ins, app_perf, color, shape in data:
    plt.plot(op_ins, app_perf / 1000, marker=shape, color=color, linestyle='None')

### Custom legend ###
# Legend elements for shapes (methods)
shape_legend = [
    Line2D([0], [0], color='black', marker='o', linestyle='None', label='Method: Give'),
    Line2D([0], [0], color='black', marker='^', linestyle='None', label='Method: Take'),
]

# Legend elements for colors (size configurations)
color_legend = [
    Line2D([0], [0], color='g', marker='o', linestyle='None', label='Size: 769 x 1024'),
    Line2D([0], [0], color='b', marker='o', linestyle='None', label='Size: 1537 x 2048'),
    Line2D([0], [0], color='orange', marker='o', linestyle='None', label='Size: 3073 x 4096'),
    Line2D([0], [0], color='r', marker='o', linestyle='None', label='Size: 6145 x 8192'),
]


# Labels and grid
plt.title(title)
plt.xlabel(x_axis_label)
plt.ylabel(y_axis_label)
plt.xticks(x_ticks, x_tick_labels)
plt.yticks(y_ticks, y_tick_labels)

# Set axis limits
plt.xlim(x_min, x_max)
plt.ylim(y_min, y_max)

# Legend
plt.legend(handles= shape_legend + color_legend, loc='lower left')

# Show plot
plt.tight_layout()
plt.show()


plt.savefig('roofline_model.png')