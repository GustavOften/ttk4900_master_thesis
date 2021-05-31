import matplotlib.pyplot as plt
import sys
import numpy as np
import re
import scipy.io


with open('sin_wave_001_5_ball_right.txt', 'r') as f:
	data = f.read()

def makere(_s):
	# return s.subs('[', '\[')
	f = r'([\deE\+\-\.]+)'
	s = _s.replace('[', R'\[')
	s = s.replace(']', R'\]')
	s = s.replace('(', R'\(')
	s = s.replace(')', R'\)')
	s = s.replace('%f', f)
	return s

m = []
k = []
l = []
expr = makere('[inf] t=%f,torque=%f,theta=%f,phi=%f,dtheta=%f,dphi=%f,x=%f,y=%f')
observer = makere('[obs] q1 : %f q2: %f dq1: %f dq2: %f')
observer_no_input = makere('[obs_no_input] q1 : %f q2: %f dq1: %f dq2: %f')

for line in data.split('\n'):
	ans = re.match(expr, line)
	obs = re.match(observer,line)
	obs_no_input = re.match(observer_no_input,line)
	if ans is None:
		if obs is None:
			if obs_no_input is None:
				continue
			q1,q2,dq1,dq2 = obs_no_input.groups()
			l += [(q1,q2,dq1,dq2)]
			continue
		q1,q2,dq1,dq2 = obs.groups()
		k += [(q1,q2,dq1,dq2)]
		continue
	t,torque,theta,phi,dtheta,dphi,x,y = ans.groups()
	m += [(t,torque,theta,phi,dtheta,dphi,x,y)]

m = np.array(m, dtype=float)
#k = np.array(k,dtype=float)
#l = np.array(l,dtype=float)
plt.plot(m[:,0], m[:,1]*100, label='torque')

plt.plot(m[:,0], m[:,2], label='theta')
plt.plot(m[:,0], m[:,3], label='phi')

plt.plot(m[:,0], m[:,4], label='dtheta')
plt.plot(m[:,0], m[:,5], label='dphi')

plt.grid()
plt.legend()
plt.show()
np.savetxt('sin_wave_001_5_ball_right_matlab.txt',m)
#np.savetxt('observer.txt',k)
#np.savetxt('observer_no_tau.txt',l)