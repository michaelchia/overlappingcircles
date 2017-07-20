#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 23:25:21 2017

@author: m
"""
from math import sqrt, atan2
from numpy import mean


def dist(coord1, coord2):
    return sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2)
    
class CircleOverlap:
    def __init__(self, circles):
        self.circles = circles
        self.vertices = None
    
    def get_vertices(self, reset = False):   
        if reset:
            self.vertices = None
        if self.vertices is None:
            self.circles.sort()
            for c in self.circles:
                self.add_circle(c)
                if self.vertices is False:
                    return False
        return self.vertices
                
    def add_circle(self, circle):
        if self.vertices is None: # first circle
            self.vertices = [circle.to_point()]
            return
        if len(self.vertices) == 1: # only a circle
            intersection = self.vertices[0].circles[0].intersect(circle)
            if intersection is None:
                self.vertices = False
                return
            if type(intersection) is Circle:
                self.vertices = [intersection.to_point()]
                return
            self.vertices = intersection
            return
        circles = self.filter_vertices(circle)
        if len(circles) == 0: # circle within vertices
            circles = self.get_circles()
            self.vertices = []
        for c in circles:
            candidates = circle.intersect(c)
            if candidates is None:
                self.vertices = False
                return
            if type(candidates) is list:
                points = self.filter_points(candidates, circles)
                self.vertices = self.vertices + points
        if len(self.vertices) == 0: # circle inside all previous circles
            self.vertices = [circle.to_point()]
        return
    
    def filter_vertices(self, circle):
        vertices = []
        circles = set()
        for pt in self.vertices:
            if circle.contains(pt.coord):
                vertices.append(pt)
                circles.add(pt.circles[0])
                circles.add(pt.circles[1])
        if len(vertices) == 0:
            return []
        self.vertices = vertices
        return circles
    
    def get_circles(self):
        circles = set()
        for pt in self.vertices:
            circles.add(pt.circles[0])
            circles.add(pt.circles[1])
        return circles
    
    def filter_points(self, candidates, circles):
        points = []
        for pt in candidates:
            switch = True
            for c in circles:
                if c not in pt.circles:
                    if not(c.contains(pt.coord)):
                        switch = False
                        break
            if switch:
                points.append(pt)
        return points
        
    # slow simple version
    def _get_vertices(self, reset = False):
        if reset:
            self.vertices = None
        if self.vertices is not None:
            return self.vertices
        # get all intersections
        candidates = []
        for i in range(0,len(self.circles)):
            c1 = self.circles[i]
            for j in range(i+1,len(self.circles)):
                c2 = self.circles[j]
                candidate = c1.intersect(c2)
                if candidate is None: # no common overlap
                    self.vertices = False
                    return False
                if type(candidate) is not Circle:
                    candidates = candidates + candidate
        if len(candidates) == 0: # a circle is contained within all other circles
            vertices = [min(self.circles).to_point()]
            self.vertices = vertices
            return vertices
        # check for intersections contained within all circles
        vertices = []
        for pt in candidates:
            switch = True
            for c in self.circles:
                if c not in pt.circles and not(c.contains(pt.coord)):
                    switch = False
                    break
            if switch:
                vertices.append(pt)
        self.vertices = vertices
        return vertices
    
    # use avg for now instead of centre of gravity
    def centroid(self):
        vertices = self.get_vertices()
        if vertices is False:
            return None
        return (mean(list(map(lambda x: x.coord[0], vertices))), \
                mean(list(map(lambda x: x.coord[1], vertices))))
    
    def order_vertices(self):
        centre = self.average()
        self.vertices.sort(key = lambda x: x.angle(centre))
        return self.vertices

    def average(self):
        vertices = self.get_vertices()
        if vertices is False:
            return None
        return (mean(list(map(lambda x: x.coord[0], vertices))), \
                mean(list(map(lambda x: x.coord[1], vertices))))
        
class Circle:
    def __init__(self, coord, r):
        self.coord = coord
        self.r = float(r)
       
    def intersect(self, other):
        d = dist(self.coord, other.coord)
        if d >= self.r + other.r or (d == 0 and self.r == other.r):
            return None
        if d <= abs(self.r - other.r):
            return min(self, other)
        l = (self.r**2 - other.r**2 + d**2) / (2*d)
        h = sqrt(self.r**2 - l**2)
        xmid = self.coord[0] + l/d*(other.coord[0] - self.coord[0])
        ymid = self.coord[1] + l/d*(other.coord[1] - self.coord[1])
        x1 = xmid + h/d*(other.coord[1] - self.coord[1])
        y1 = ymid - h/d*(other.coord[0] - self.coord[0])
        x2 = xmid - h/d*(other.coord[1] - self.coord[1])
        y2 = ymid + h/d*(other.coord[0] - self.coord[0])
        return [IntersectionPoint((x1,y1),(self,other)), \
                IntersectionPoint((x2,y2),(self,other))]
        
    def contains(self, coord):
        return self.r > dist(self.coord, coord)
    
    def to_point(self):
        return IntersectionPoint(self.coord, [self])
    
    def __hash__(self):
        return hash(hash(self.coord) + hash(self.r))
    def __eq__(self, other):
        return self.coord == other.coord and self.r == other.r
    def __gt__(self, other):
        return self.r > other.r

class IntersectionPoint:
    def __init__(self, coord, circles):
        self.coord = coord
        self.circles = circles
    def angle(self, centre):
        return atan2(self.coord[1]-centre[1], self.coord[0]-centre[0])
    def __hash__(self):
        return hash(self.coord)
    def __eq__(self, other):
        eps = 1e-10
        x_diff = abs(self.coord[0] - other.coord[0])
        y_diff = abs(self.coord[1] - other.coord[1])
        return x_diff <= eps and y_diff <= eps
    def __str__(self):
        return str(self.coord)
    def __repr__(self):
        return self.__str__()
