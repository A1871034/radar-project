from PIL import Image

if __name__ == "__main__":
    im = Image.open(r"Figures/Wilpena_Pound_Cross_Section.jpg")
    data = im.getdata()
    
    res = b""

    for n in range(im.size[0]):
        for m in range(im.size[1]):
            if data[m*im.size[0] + n] != 255:
                res += (im.size[1]-m-1).to_bytes(4, 'little')
                break

    print(f"{len(res)//4}/{im.size[0]}")

    with open("2D/heights.data", "wb") as f:
        f.write(res)
    
