import matplotlib.pyplot as plt
import matplotlib.animation as animation
import matplotlib.colors
import matplotlib.image as image
import numpy as np

#from IPython import display 

if __name__ == "__main__":
    with open("2D/sim.data", "r") as f:
        header = f.readline().strip().split(',')
        seekTo = f.tell()
    SIZE_X, SIZE_Y, PRECISION_BYTES = [int(field) for field in header]

    print(header)
        
    with open("2D/sim.data", "rb") as f:
        f.seek(seekTo)
        data = np.frombuffer(f.read(), dtype=np.float)
        finalTell = f.tell()

    STEPS = (finalTell-seekTo) // (SIZE_X*SIZE_Y*PRECISION_BYTES)
    #V_MIN, V_MAX = np.min(data), np.max(data)
    #v_max = max(abs(V_MIN), V_MAX)

    data = data.reshape((STEPS, SIZE_X, SIZE_Y)).swapaxes(1,2)

    v_max = 0.06

    terrain = image.imread(r"Figures/Wilpena_Pound_Cross_Section.jpg")

    fig = plt.figure(figsize=(SIZE_X//100,SIZE_Y//100))
    im = plt.imshow(data[0], interpolation='none', cmap="coolwarm", aspect=1, vmin=-v_max, vmax=v_max)
    ax = plt.gca()
    ax.invert_yaxis()
    #ax.imshow(terrain)
    ax.set(xlabel="index (n)", ylabel="index (m)")
    #cbar = plt.colorbar(im)
    ax.set_title("Electric Field Strength at Step = 0")

    # Uncomment below for running avg scaling
    m_val = [v_max]
    def update(t):
        im.set_array(data[t])

        # Code below will adjust scale for better dynamic range while 
        # maintaining relative data between positive and negative values

        ax.set_title(f"Electric Field Strength at Step = {t}")
        return [im]
    
    FPS = 30
    ani = animation.FuncAnimation(fig=fig, func=update, frames=STEPS, interval=100/FPS)

    print("Saving animation")
    
    plt.show()

    #ani.save("tmp/2D_FDTD_Example.gif", "ffmpeg", FPS)

    