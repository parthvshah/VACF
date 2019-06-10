f = open('HISTORY', 'r')
lines = f.readlines()
for x in range(len(lines)):
    if(lines[x][0]=='t'):
        time = int(lines[x][8:25])
    if(lines[x][0]=='A'):
        particle = int(lines[x][2:20])
        print(str(time)+" "+str(particle)+" "+lines[x+2].strip())