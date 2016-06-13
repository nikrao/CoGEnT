from numpy import *
import scipy as Sci
import scipy.linalg

def find_next_atom_l1(gradf):
    s = sign(-gradf);
    p = gradf.shape;
    p = p[0];
    idx = argmax(absolute(gradf));
    atom = zeros((p,1));
    atom[idx] = s[idx];
    return(atom);