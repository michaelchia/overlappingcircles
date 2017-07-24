#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jul 21 00:46:13 2017

@author: m
"""
from math import pi, sin, cos, sqrt
from random import uniform, randint
import matplotlib.pyplot as plt
from circleoverlap import CircleOverlap, Circle, dist, find_overlaps, \
    max_overlap

# Test functions
def generate_overlap_circles(n, 
                              centre = (0,0), 
                              rng = (1,10), 
                              min_overlap = 0.1, 
                              r_transform = lambda x: x):
    x, y = centre
    circles = []
    for i in range(0,n):
        r = r_transform(uniform(*rng))
        d = uniform(0,r - min_overlap)
        theta = uniform(0,2*pi)
        circles.append(Circle((d*cos(theta)+x,d*sin(theta)+y), r)) 
    return circles

def generate_circles(n, xy_rng = (-10, 10), r_rng = (1,10), r_transform = lambda x: x):
    circles = []
    for i in range(0,n):
        x = uniform(*xy_rng)
        y = uniform(*xy_rng)
        r = r_transform(uniform(*r_rng))
        circles.append(Circle((x,y), r))
    return circles

def plot_circles(circles_ls, 
                 rng = (-20, 20), 
                 circles = True,
                 circle_points = False,
                 overlaps = True,
                 vertices = True,
                 centroid = True,
                 component_centroids = False):
    fig, ax = plt.subplots()
    #ax.set_xlim(rng)
    #ax.set_ylim(rng)
    if circles:
        for c in circles_ls:
            ax.add_artist(plt.Circle(c.coord(), c.r, color='0.2', fill=False))
    if overlaps:
        overlap_ls = find_overlaps(circles_ls)
        for overlap in overlap_ls:
            plot_overlap(overlap, subplots = (fig, ax), 
                         circles = False,
                         circle_points = circle_points,
                         vertices = vertices,
                         centroid = centroid,
                         component_centroids = component_centroids)
    
def plot_overlap(overlap, 
                 subplots = None,
                 rng = (-20, 20),
                 circles = True,
                 circle_points = False,
                 vertices = True,
                 centroid = True,
                 component_centroids = False):
    if subplots is None:
        subplots = plt.subplots()       
    #    subplots[1].set_xlim(rng)
    #    subplots[1].set_ylim(rng)
    fig, ax = subplots
    if circles:
        for c in overlap.circles:
            ax.add_artist(plt.Circle(c.coord(), c.r, color='0.2', fill=False))
    if circle_points:
        for c in overlap.circles:
            ax.plot(c.x,c.y,'o',color='0.2',ms=2) 
    if vertices:
        for v in overlap.get_vertices():
            ax.plot(v.x,v.y,'ob',ms=3)
    if centroid:
        cent = overlap.centroid()
        ax.plot(cent[0],cent[1],'or',ms=3)
    if component_centroids:
        if len(overlap.get_vertices()) > 2:
            for seg in overlap.get_circle_segments():
                cent = seg.centroid()
                ax.plot(cent[0],cent[1],'xr')
    return (fig, ax)
        
def compare_vertices(vert1, vert2):
    eps = 1e-12
    for v1 in vert1:
        match = False
        for v2 in vert2:
            x_diff = abs(v1.x - v2.x)
            y_diff = abs(v1.y - v2.y)
            if x_diff <= eps and y_diff <= eps:
                match = True
                break
        if not(match):
            return False
    return True

def test_vertices(iterations, rng=(1,100)):
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
    return fails