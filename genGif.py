import sys
from PIL import Image, ImageDraw

images = []

input_file = sys.argv[1]
output_file = sys.argv[2]

height = 0
width = 0
step = 0
with open(input_file, 'r') as f:
    line = f.readline().strip().split()
    height = int(line[0])
    width = int(line[1])
    step = int(line[2])

    for s in range(step):
        im = Image.new('RGB', (width, height), (0,0,0))
        for i in range(height):
            line = list(map(int, f.readline().strip().split()))
            for j, bolt in enumerate(line):
                if bolt > 0:
                    l = int((bolt - 0.5) * 0.5 * 128)
                    im.putpixel((j, i), (l, l, l))
        images.append(im)
        f.readline()

if len(images) > 1:
    images[0].save(output_file,
               save_all=True, append_images=images[1:], optimize=False, duration=500, loop=0)