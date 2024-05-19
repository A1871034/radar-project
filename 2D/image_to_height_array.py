from PIL import Image

if __name__ == "__main__":
    im = Image.open(r"Figures/Wilpena_Pound_Cross_Section.jpg")
    data = im.getdata()
    
    res = b""
    for n in range(im.size[0]):
        first = None
        last = im.size[1] - 1
        for m in range(im.size[1]):
            if first is None:
                if data[m*im.size[0] + n] != 255:
                    first = (im.size[1]-m-1)
            else:
                if data[m*im.size[0] + n] == 255:
                    last = (im.size[1]-m-1)

        res += ((first+last)//2).to_bytes(4, 'little')
				

    print(f"{len(res)//4}/{im.size[0]}")

    with open("2D/heights.data", "wb") as f:
        f.write(res)
    
