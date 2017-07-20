#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 00:46:13 2017

@author: m
"""

from random import uniform, randint
import matplotlib.pyplot as plt
import time
from circleoverlap import CircleOverlap, Circle, dist

# Testing
def generate_circles(n, rng = (-10, 10), centre = (0,0), overlap_max = 5):
    circles = []
    for i in range(0,n):
        x = uniform(rng[0],rng[1])
        y = uniform(rng[0],rng[1])
        r = dist((x,y), centre) + uniform(0, overlap_max)
        circles.append(Circle((x,y), r))
    return circles

def plot_circles(circles, rng = (-20, 20)):
    fig, ax = plt.subplots()
    ax.set_xlim(rng)
    ax.set_ylim(rng)
    for c in circles:
        ax.add_artist(plt.Circle(c.coord, c.r, color='b', fill=False))
        
def compare_vertices(vert1, vert2):
    for v1 in vert1:
        match = False
        for v2 in vert2:
            if v1 == v2:
                match = True
                break
        if not(match):
            return False
    return True

# test vertice algo
iterations = 1000
count = 0
for i in range(0,iterations):   
    circles = generate_circles(randint(1,100))
    overlap = CircleOverlap(circles)
    vert1 = overlap.get_vertices(reset = True)
    vert2 = overlap._get_vertices(reset = True)
    if compare_vertices(vert1,vert2):
        count = count + 1
    print(i+1)
print(count == iterations)

# time vertice algo
circles = generate_circles(4)
overlap = CircleOverlap(circles)
plot_circles(circles)
t0 = time.time()
vert = overlap.get_vertices(reset = True)
t1 = time.time()
print(t1-t0)
print(vert)