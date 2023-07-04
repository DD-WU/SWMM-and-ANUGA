
x=[i for i in range(0,80,1)]
with open("test.txt","w") as f:
    f.write(str(x))
f.close()