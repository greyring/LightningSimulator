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
    power = 5
    lightning_point = width-1
lightning_width = int(width/(lightning_point+1))

with open(output_file, 'w') as f:
    f.write(str(width)+" "+str(width)+" "+str(power)+" "+str(eta)+"\n")
    f.write(str(lightning_point)+"\n")
    if mode == 1:
        f.write(str(mid_point)+" "+str(mid_point)+"\n")
        total_num = (width-1)*4
        f.write(str(total_num)+"\n")
        for i in range(width-1):
            f.write("0 "+str(i)+"\n")
            f.write(str(i+1)+" 0\n")
            f.write(str(width-1)+" "+str(i+1)+"\n")
            f.write(str(i)+" "+str(width-1)+"\n")
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
        ground_num = 10
        grount_interval = int(width/(ground_num+1))
        ground_col = grount_interval
        f.write(str(ground_num)+"\n");
        for i in range(1, ground_num+1):
            height = int(random.uniform(0, 1)*width/5+width*4/5)
            height = min(height, width-1)
            cur_col = i*grount_interval
            f.write(str(height)+" "+str(cur_col)+"\n")
            
    
