import numpy as np
import matplotlib.pyplot as plt

n = [100, 200, 300, 500, 1000, 5000, 10000, 20000]
julia_times =  [26.7, 43.6, 64.1, 106.2, 200, 987, 1.96*1000, 3.91*1000]
python_times = [142, 239, 368, 607, 1.17*1000, 5.86*1000, 11.6*1000, 23.5*1000]

N = len(n)

fig, ax = plt.subplots()

ind = np.arange(N)    # the x locations for the groups
width = 0.35         # the width of the bars

p1 = ax.bar(ind, julia_times, width, bottom=0)
p2 = ax.bar(ind + width, python_times, width, bottom=0)

ax.set_xticks(ind + width / 2)
ax.set_xticklabels(n)
ax.set_xlabel('Number of nodes')

ax.legend((p1[0], p2[0]), ('Julia', 'Python'))
ax.set_yscale('log')
ax.autoscale_view()
ax.set_ylabel('Runtime/μs')

plt.show()