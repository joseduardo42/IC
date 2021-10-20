from math import sin, cos, pi
import numpy as np
from numpy.linalg import inv

from Ex4_harmonicbalance_nonlin_circ import F_inv

h = 2

def F(f):
  F = np.array([[1, sin(2*pi*0*1/(2*h+1)), cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2)/(2*h+1))],
                [1, sin(2*pi*1*1/(2*h+1)), cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2)/(2*h+1))],
                [1, sin(2*pi*2*1/(2*h+1)), cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2)/(2*h+1))],
                [1, sin(2*pi*3*1/(2*h+1)), cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2)/(2*h+1))],
                [1, sin(2*pi*4*1/(2*h+1)), cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2)/(2*h+1))]])

  F_adiantada = F + 2*pi*f
  return[F_adiantada]

def F_inv():
  F = np.array([[1, sin(2*pi*0*1/(2*h+1)), cos((2*pi*0*1/(2*h+1))), sin(2*pi*0*2/(2*h+1)), cos((2*pi*0*2)/(2*h+1))],
                [1, sin(2*pi*1*1/(2*h+1)), cos((2*pi*1*1/(2*h+1))), sin(2*pi*1*2/(2*h+1)), cos((2*pi*1*2)/(2*h+1))],
                [1, sin(2*pi*2*1/(2*h+1)), cos((2*pi*2*1/(2*h+1))), sin(2*pi*2*2/(2*h+1)), cos((2*pi*2*2)/(2*h+1))],
                [1, sin(2*pi*3*1/(2*h+1)), cos((2*pi*3*1/(2*h+1))), sin(2*pi*3*2/(2*h+1)), cos((2*pi*3*2)/(2*h+1))],
                [1, sin(2*pi*4*1/(2*h+1)), cos((2*pi*4*1/(2*h+1))), sin(2*pi*4*2/(2*h+1)), cos((2*pi*4*2)/(2*h+1))]])

  inverse_F = inv(F)
  return[inverse_F]

def omega(w):
  #omega matrix
  omega = np.zeros((5,5))
  omega[1, 2] = -w; omega[2, 1] = w; omega[3, 4] = -2*w; omega[4, 3] = 2*w
  return[omega]
