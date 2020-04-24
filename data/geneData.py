import sys
import random 

output_file = sys.argv[1]
width = int(sys.argv[2])
mode = int(sys.argv[3])
power = 1
eta = 1
mid_point = int(width/2)
if mode == 1:
    lightning_point = 1
elif mode == 2:
    lightning_point = 1
else:
    lightning_point = 5
lightning_width = int(width/(lightning_point+1))

with open(output_file, 'w') as f:
    f.write(str(width)+" "+str(width)+" "+str(power)+" "+str(eta)+"\n")
    f.write(str(lightning_point)+"\n")
    if mode == 1:
        f.write("0 "+str(mid_point)+"\n")
        f.write("1\n");
        f.write(str(width)+" "+str(mid_point)+"\n")
    elif mode == 2:

        f.write("0 "+str(mid_point)+"\n")
        f.write(str(width)+"\n");
        for i in range(width):
            height = int(random.uniform(0, 1)*width/4+width*3/4)
            f.write(str(height-1)+" "+str(i)+"\n")
    elif mode == 3:
        lightning_col = lightning_width
        for i in range(lightning_point):
            f.write("0 "+str(lightning_col)+"\n")
            lightning_col += lightning_width
        f.write(str(width)+"\n");
        for i in range(width):
            height = int(random.uniform(0, 1)*width/5+width*4/5)
            height = max(height, width-1)
            f.write(str(height)+" "+str(i)+"\n")
    
