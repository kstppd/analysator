#!/usr/bin/python
# 
# This file is part of Analysator.
# Copyright 2013-2016 Finnish Meteorological Institute
# Copyright 2017-2018 University of Helsinki
# 
# For details of usage, see the COPYING file and read the "Rules of the Road"
# at http://www.physics.helsinki.fi/vlasiator/
# 
# This program is free software; you can redistribute it and/or modify
# it under the terms of the GNU General Public License as published by
# the Free Software Foundation; either version 2 of the License, or
# (at your option) any later version.
# 
# This program is distributed in the hope that it will be useful,
# but WITHOUT ANY WARRANTY; without even the implied warranty of
# MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
# GNU General Public License for more details.
# 
# You should have received a copy of the GNU General Public License along
# with this program; if not, write to the Free Software Foundation, Inc.,
# 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA.
# 

import pytools as pt
import numpy as np
import sys
from multiprocessing import Pool
from tqdm import tqdm
RE=6378137.0

def mag(v):
    return np.sqrt(np.dot(v,v))

def unit(v):
    return v/mag(v)

def cash_karp(origin,f,data,var,ds,tol):

    #k1
    b=f.read_interpolated_fsgrid_variable(var,origin, operator="pass",periodic=[True,True,True])
    k1=unit(b)

    #k2
    step=origin+k1*ds/5.0;
    b=f.read_interpolated_fsgrid_variable(var,step, operator="pass",periodic=[True,True,True])
    k2=unit(b)

    #k3
    step=origin+ k1*(3.0*ds/40.0)+k2*(9.0*ds/40.0)
    b=f.read_interpolated_fsgrid_variable(var,step, operator="pass",periodic=[True,True,True])
    k3=unit(b)

    #k4
    step=origin+k1*(3.0*ds/10.0)+k2*(-9.0*ds/10.0)+k3*(6.0*ds/5.0)
    b=f.read_interpolated_fsgrid_variable(var,step, operator="pass",periodic=[True,True,True])
    k4=unit(b)

    #k5
    step=origin+k1*(-11.0*ds/54.0)+k2*(5.0*ds/3.0)+k3*(-70.0*ds/27.0)+k4*(35.0*ds/27.0)
    b=f.read_interpolated_fsgrid_variable(var,step, operator="pass",periodic=[True,True,True])
    k5=unit(b)

    #k6
    step=origin+k1*(1631.0*ds/55296.0)+k2*(175.0*ds/512.0)+k3*(575.0*ds/13824.0)+k4*(44275.0*ds/110592.0)+k5*(253.0*ds/4096)
    b=f.read_interpolated_fsgrid_variable(var,step, operator="pass",periodic=[True,True,True])
    k6=unit(b)

    #estimates
    rk4=np.zeros(3)
    rk5=np.zeros(3)
    rk4=origin+ k1*2825.0*ds/27648.0 +k3*18575.0*ds/48384.0+ k4*13525.0*ds/55296.0+ k5*277.0*ds/14336.0+ k6*1.0*ds/4.0;
    rk5=origin+ k1*37.0*ds/378.0 +k3*250.0*ds/621.0+ k4*125.0*ds/594.0+ k6*512.0*ds/1771.0;
  
    #error & step handling
    err= 100.0*abs((rk5-rk4)/rk4)
    maxErr=np.max(err)
    S =0.95*np.power(tol/maxErr, 0.2);
    ds=S*ds
    if (maxErr>tol):
        return False;

    origin=rk5
    return True

def sphere(r):
    R = r*RE
    tht = np.linspace(0, 2* np.pi, 100)
    phi = np.linspace(0, 2* np.pi, 100)
    tht, phi = np.meshgrid(tht, phi)
    x = (R * np.cos(tht))* np.cos(phi)
    y = (R * np.cos(tht))* np.sin(phi)
    z = R * np.sin(tht)
    return np.array([x, y, z])

def curl_on_sphere(file,r):
    f=pt.vlsvfile.VlsvReader(file)
    b=f.read_fsgrid_variable("fg_b")
    ds=1e4
    curls=[]

    for r0 in tqdm(r):
        x,y,z=r0

        #x-component (Fz/theta_y -Fy/theta_z)
        rp=[x,y+ds,z]
        rm=[x,y-ds,z]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[2]-bm[2])/(2*ds);

        rp=[x,y,z+ds]
        rm=[x,y,z-ds]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[1]-bm[1])/(2*ds);
        curlBx=term1-term2;

        #y-component (Fx/theta_z -Fz/theta_x)
        rp=[x,y,z+ds]
        rm=[x,y,z-ds]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[0]-bm[0])/(2*ds);

        rp=[x+ds,y,z]
        rm=[x-ds,y,z]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[2]-bm[2])/(2*ds);
        curlBy=term1-term2;

        #z-component (Fy/theta_x -Fx/theta_y)
        rp=[x+ds,y,z]
        rm=[x-ds,y,z]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[1]-bm[1])/(2*ds);

        rp=[x,y+ds,z]
        rm=[x,y-ds,z]
        bp=f.read_interpolated_fsgrid_variable("fg_b",rp, operator="pass",periodic=[True,True,True])
        bm=f.read_interpolated_fsgrid_variable("fg_b",rm, operator="pass",periodic=[True,True,True])
        term1= (bp[0]-bm[0])/(2*ds);
        curlBz=term1-term2;

        curl=np.array([curlBx,curlBy,curlBz])
        curl_mag=mag(curl)
        curls.append(curl_mag)
    return curls

def detect_cusps(file):
    r=sphere(7)
    curls=curl_on_sphere(file,r)


def main():
    files=sys.argv[1::]
    for file in files:
        detect_cusps(file)



main()



