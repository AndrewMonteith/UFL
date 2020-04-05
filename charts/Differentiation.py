import numpy as np
import matplotlib.pyplot as plt

n = [1, 2, 3, 4, 5]
julia_times =  [458.7/1000, 1.82, 1.85, 2.96, 2.98]
python_times = [1.21, 3.51, 3.45, 5.05, 5.74]

N = len(n)

fig, ax = plt.subplots()

ind = np.arange(N)    # the x locations for the groups
width = 0.35         # the width of the bars

p1 = ax.bar(ind, julia_times, width, bottom=0)
p2 = ax.bar(ind + width, python_times, width, bottom=0)

ax.set_xticks(ind + width / 2)
ax.set_xticklabels(n)
ax.set_xlabel('Instance')

# ax.set_yscale('log')
ax.autoscale_view()
ax.set_ylabel('Runtime/ms')

ax.legend((p1[0], p2[0]), ('Julia', 'Python'))

plt.show()