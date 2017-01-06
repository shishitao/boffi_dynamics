# ad hoc script to compute the funny colormap used by "map.gp"

values = ( 1.25,     2.0,    3.0,      4.0,    5.0 )
colors = ('yellow', 'cyan', 'yellow', 'cyan', 'red')

begins = (255, 255, 255)
ends   = (245, 150, 155)
maplen = 10.0
delta  = map(lambda x, y: (y-x)/maplen, begins, ends)

def color(x):
    return tuple(map(int, map(lambda s,d: s+d*x, begins, delta)))

print """\
set cbrange [0:10]
set cntrparam levels discrete""",
for v in values: print "%f,"%(v,),
print 10.0
print 'set palette defined ( 0.000 "#%2.2x%2.2x%2.2x",' % begins,

for v, c in zip(values, colors):
    x1, x2 = v*0.95, v*1.05
    c1, c2 = color(x1), color(x2)
    print "%4.3f" % x1,
    print '"#%2.2x%2.2x%2.2x",' % c1,
    print '%4.3f "%s", ' % (v, c),
    print "%4.3f" % x2,
    print '"#%2.2x%2.2x%2.2x",' % c2,

print '10.000 "#%2.2x%2.2x%2.2x")' % ends
