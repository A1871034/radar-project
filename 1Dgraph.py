from matplotlib import pyplot as plt
import matplotlib.animation as animation
import pandas as pd
import numpy as np

df = pd.read_csv("test.csv", header=None)

SAVE = False

print("Animating")

fig = plt.figure(figsize=(7,7))
ax = plt.axes(xlim=(0,199),ylim=(-1,1))
scat=ax.plot(df.columns, list(df.iloc[0].to_numpy()), label="t=0 (steps)")[0]
ax.set(xlabel="index (i)", ylabel="Ez[i] (V\m)", title="1D FDTD Demonstration", aspect=75)
ax.legend()


def update(t):
    scat.set_xdata(df.columns)
    scat.set_ydata((df.iloc[t].to_numpy()))
    scat.set_label(f"t={t*3} (steps)")
    ax.legend()
    #scat.set_offsets([[i, df[i][t]] for i in range(200)])
    return scat

ani = animation.FuncAnimation(fig=fig, func=update, frames=len(df), interval=50)
print("Displaying")
plt.show()

if SAVE:
    print("Saving")
    ani.save(filename="tmp/1d_FDTD_example.gif", writer="pillow")
