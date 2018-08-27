import numpy as np



pi=3.141592653586793238462643383279502884197









#a=np.array([1,2,3])
#b=np.array([4,5,6])

#def f(x):
#	return np.roll(a,x,0)

#t=np.array([f(0).dot(b),f(1).dot(b),f(2).dot(b)])

#print np.argmax(t)
















a=np.array([0,.5,1,.5,0,.5,1,.5,0])
b=np.array([1,.5,0,.5,1,.5,0,.5,1])


def f(x):
	return np.roll(a,x,0)


def g(s,t):
	k = []
	for i in xrange(len(s)):
		k.append(f(i).dot(b))
	return k

c = g(a, b)
print c

print np.argmax(c)









