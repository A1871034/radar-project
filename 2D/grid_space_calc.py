from numpy import ceil

if __name__ == "__main__":
    unit = input("unit: ").lower()
    freq = int(input("freq: "))
    ppw = int(input("grid points per wavelength: "))

    unit_map = {
        "ghz": 10**9,
        "mhz": 10**6,
        "hz": 1
    }

    if unit.lower() not in unit_map:
        print("Unkown unit")
        exit()

    freq = unit_map[unit] * freq
    c = 3 * 10**8

    dx_and_dy = c/freq/ppw

    print(f"Δx = Δy = {dx_and_dy} meters")
    print(f"Δt = {1/(freq*ppw)} seconds") # freq*ppw = dx_and_dy/c

    repeat = "y"
    while repeat == "y":
        des_width = int(input("Desired width (m): "))
        des_height = int(input("Desired height (m): "))

        width = int(ceil(des_width/dx_and_dy))
        height = int(ceil(des_height/dx_and_dy))

        print(f"Width (px) = {width}")
        print(f"Height (px) = {height}")
        
        precision = 8 # bytes
        print(f"Assuming {precision}B precision")


        mem_usage = 0 # Bytes
        
        mem_usage += width * (height - 1) * precision # hx
        mem_usage += (width -1 ) * height * precision # hy
        mem_usage += width * height * precision       # ez
        mem_usage += width * 6 * 2 * precision        # top + bot ABC buffers
        mem_usage += height * 6 * 2 * precision       # left + right ABC buffers

        print(f"Grid Memory Usage: {mem_usage/(10**9)} GB")
        repeat = input("Re-input desired dimensions? (Y/N): ").lower()
