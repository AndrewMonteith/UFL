import numpy as np
import matplotlib.pyplot as plt

n = [100, 200, 300, 500, 1000, 5000, 10000, 20000]
julia_times_pre_order =  [12.3, 22.1, 34.6, 54.3, 107, 522.6, 1*1000, 2.1*1000]
python_times_pre_order = [19.8, 40.4, 63.4, 112, 233, 1.18*1000, 2.42*1000, 4.87*1000]

julia_times_post_order = [18.8, 35.4, 50.8, 88.2, 171, 820, 1.65*1000, 3.26*1000]
python_times_post_order = [79.7, 161, 250, 449, 915, 4.75*1000, 9.61*1000, 20.3*1000]

N = len(n)

fig, (ax, ax2) = plt.subplots(ncols=2)

ind = np.arange(N)    # the x locations for the groups
width = 0.35         # the width of the bars

ax.set_title("Pre-Order Traversal")
ax2.set_title("Post-Order Traversal")

p1 = ax.bar(ind, julia_times_pre_order, width, bottom=0)
p2 = ax.bar(ind + width, python_times_pre_order, width, bottom=0)

p3 = ax2.bar(ind, julia_times_post_order, width, bottom=0)
p4 = ax2.bar(ind + width, python_times_post_order, width, bottom=0)

for chart in [ax, ax2]:
    chart.set_xticks(ind + width / 2)
    chart.set_xticklabels(n)
    chart.set_xlabel('Number of nodes')

    chart.set_yscale('log')
    chart.autoscale_view()
    chart.set_ylabel('Runtime/μs')

ax.legend((p1[0], p2[0]), ('Julia', 'Python'))
ax2.legend((p3[0], p4[0]), ('Julia', 'Python'))

plt.show()