data = 'C = (39.447, 94.657, 11.824) N = (39.292, 95.716, 11.027) Ca = (39.462, 97.101, 11.465)'
t1 = data.replace ('(', ',') # change ( to ,
t2 = t1.replace (')', ',')    # change ) to ,
l  = t2.split (',')           # split on ,

p = (float(l[1]),float(l[2]), float(l[3]))
q = (float(l[5]),float(l[6]), float(l[7]))
r = (float(l[9]),float(l[10]), float(l[11]))
print(p)
print(q)
print(r)