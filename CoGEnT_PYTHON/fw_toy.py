from numpy import *
import atom_selections as asel
import matplotlib.pyplot as plt


# toy code to see if I can run a simple frank wolfe method for sparse recovery

n = 2048; # number of observations
p = 5000; # dimension

# set random number generator for repeatibility
random.seed(10);

# data matrix
X = random.normal(0,1/sqrt(n),[n,p]);
matsize = X.shape;
print "The size of the data is %d x %d \n" %(matsize[0], matsize[1]);

# target sparse vector
wtrue = zeros((p,1));
inds  = random.randint(0,p+1,[100,1]);
isize = inds.shape;
print "the size of inds is %d %d \n" %(isize[0],isize[1]);

wtrue[inds] = 1;
print "norm of wtrue is %f \n" %(linalg.norm(wtrue));

# observations
y = dot(X,wtrue);


# FW parameters
maxiter = 1000;
tau = sum(absolute(wtrue));
w = zeros((p,1));
w[1] = tau;


# FW METHOD
objective = zeros((maxiter-1,1))
for t in range(0,maxiter-1):
      # select atom
      resid = dot(X,w)-y;
      gradf = dot(transpose(X),resid);
      atom = asel.find_next_atom_l1(gradf);
      # line search to get proper step size
      temp =y - tau*dot(X,atom);
      rw = -resid-temp;
      step = dot(transpose(-resid),rw)/dot(transpose(rw),rw);
      #step = 2/(3+t);
      w = (1-step)*w + step*atom*tau;
      
      objective[t] = linalg.norm(y - dot(X,w)) * linalg.norm(y - dot(X,w)) * 0.5;
      
      print "Objective function value %f \n" %(objective[t]);
      
# plot the objective using matplotlib
#plt.plot(objective);
#plt.ylabel('objective function value');
#plt.xlabel('iterations')
#plt.show();
      
