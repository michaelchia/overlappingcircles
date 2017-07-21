#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 00:46:13 2017

@author: m
"""
from random import uniform, randint
import matplotlib.pyplot as plt
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

def plot_circles(overlap, rng = (-20, 20)):
    fig, ax = plt.subplots()
    ax.set_xlim(rng)
    ax.set_ylim(rng)
    for c in overlap.circles:
        ax.add_artist(plt.Circle(c.coord, c.r, color='b', fill=False))
    for v in overlap.get_vertices():
        ax.plot(v.coord[0],v.coord[1],'ob')
    cent = overlap.centroid()
    #avg = overlap.average()
    #ax.plot(avg[0],avg[1],'xg')
    ax.plot(cent[0],cent[1],'or')
    if len(overlap.get_vertices()) > 2:
        for seg in overlap.get_circle_segments():
            cent = seg.centroid()
            ax.plot(cent[0],cent[1],'xr')
    #cent = Polygon(overlap.get_vertices()).centroid()
    #ax.plot(cent[0],cent[1], 'xg')
        
def compare_vertices(vert1, vert2):
    eps = 1e-12
    for v1 in vert1:
        match = False
        for v2 in vert2:
            x_diff = abs(v1.coord[0] - v2.coord[0])
            y_diff = abs(v1.coord[1] - v2.coord[1])
            if x_diff <= eps and y_diff <= eps:
                match = True
                break
        if not(match):
            return False
    return True

# test vertice algo
iterations = 1000
count = 0
fails = []
for i in range(0,iterations):   
    circles = generate_circles(randint(1,100))
    overlap = CircleOverlap(circles)
    vert1 = overlap.get_vertices(reset = True)
    vert2 = overlap._get_vertices(reset = True)
    if not(compare_vertices(vert1,vert2)):
        fails.append(overlap)
    print(i+1)
print(len(fails) == 0)

# time vertice algo
circles = generate_circles(10)
overlap = CircleOverlap(circles)
plot_circles(overlap)
print(overlap.centroid())