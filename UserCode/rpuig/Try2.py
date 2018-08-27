import numpy as np

#pi=3.14159265358979323846264338327950841971

a=np.array([1,10,4,7,2,3,9,8,5])
b=np.array([6,2,3])


def f(x):
	return np.roll(np.absolute(a),x,0)

def g(s,t):
	k = []
	for i in xrange(len(t)):
		k.append(s[i]*b[i])
	return k

c=g(np.absolute(a),np.absolute(b))


def h(y,z):
	l=[]
	for d in xrange(len(y)-len(z)+1):
		l.append(np.sum(g(f(0-d),b)))
	return l

r=h(np.absolute(a),np.absolute(b))

print r


#print max(r)

print np.argmax(r)


