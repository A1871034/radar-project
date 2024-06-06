import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors
import matplotlib.image as image
import numpy as np

if __name__ == "__main__":
    PRECISION_BYTES = 8        
    with open("sim.data", "rb") as f:
        data = np.frombuffer(f.read(), dtype=np.float)
        finalTell = f.tell()
    STEPS = (finalTell) // (PRECISION_BYTES)
    print(f"loaded {STEPS} steps")

    y_ABS_bound = max(max(data), -min(data))*1.1

    SAVE = True

    SAMPLE_RATE = 8

    print("Animating")

    fig = plt.figure(figsize=(7,7))
    ax = plt.axes(xlim=(0,STEPS*SAMPLE_RATE),ylim=(-y_ABS_bound,y_ABS_bound))
    scat=ax.plot([i*SAMPLE_RATE for i in range(STEPS)], data, label="ez[100][150])")[0]
    ax.set(xlabel="Step", ylabel="Ez[i] (V\m)", title="Electric Field Strength", aspect=STEPS/y_ABS_bound)
    ax.legend()

    if SAVE:
        print("Saving")
        plt.savefig('./tmp/Reciever.png', dpi=600)

    print("Displaying")
    plt.show()

    