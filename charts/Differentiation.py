import numpy as np
import matplotlib.pyplot as plt

n = [100, 200, 300, 500, 1000, 5000, 10000, 20000]
julia_times =  [103.5, 181, 259, 403, 761, 3.5*1000, 7*1000, 9*1000]
python_times = [348, 512, 710, 1.11*1000, 2070, 9890, 20000, 39600]

N = len(n)

fig, ax = plt.subplots()

ind = np.arange(N)    # the x locations for the groups
width = 0.35         # the width of the bars

p1 = ax.bar(ind, julia_times, width, bottom=0)
p2 = ax.bar(ind + width, python_times, width, bottom=0)

ax.set_xticks(ind + width / 2)
ax.set_xticklabels(n)
ax.set_xlabel('Number of nodes')

ax.set_yscale('log')
ax.autoscale_view()
ax.set_ylabel('Runtime/μs')

ax.legend((p1[0], p2[0]), ('Julia', 'Python'))

plt.show()