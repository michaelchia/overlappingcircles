#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 23:25:21 2017

@author: m
"""
from math import sqrt, atan2, sin, cos, pi
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
            if circle.contains(pt):
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
                    if not(c.contains(pt)):
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
                if type(candidate) is list:
                    candidates = candidates + candidate
        if len(candidates) == 0: # all circles are contained within other circles
            vertices = [min(self.circles).to_point()]
            self.vertices = vertices
            return vertices
        # check for intersections contained within all circles
        vertices = []
        for pt in candidates:
            switch = True
            for c in self.circles:
                if c not in pt.circles and not(c.contains(pt)):
                    switch = False
                    break
            if switch:
                vertices.append(pt)
        if len(vertices) == 0: # a circles is contained all within other circles
            vertices = [min(self.circles).to_point()]
            self.vertices = vertices
            return vertices
        self.vertices = vertices
        return vertices
    
    def centroid(self, area = False):
        vertices = self.get_vertices()
        if vertices is None:
            return None
        if len(vertices) == 1: # if just a circle
            c = vertices[0].circles[0].coord
            if area:
                return (c, vertices[0].circles[0].area())
            return c
        if len(vertices) == 2: # if only two intersecting circles
            s_seg = min(vertices[0].circles).to_segment(vertices[0],vertices[1])
            l_seg = max(vertices[0].circles).to_segment(vertices[0],vertices[1])
            if l_seg.inv_contains(s_seg):
                s_c, s_a = s_seg.inv_centroid(area = True)
            else:
                s_c, s_a = s_seg.centroid(area = True)
            l_c, l_a = l_seg.centroid(area = True)
            a = s_a + l_a
            x = (s_c[0] * s_a + l_c[0] * l_a) / a
            y = (s_c[1] * s_a + l_c[1] * l_a) / a
            if area:
                return ((x, y), a)
            return (x, y)
        shapes = self.get_circle_segments()
        shapes.append(Polygon(vertices))
        a, x, y = (0, 0, 0)
        for s in shapes:
            cent, ar = s.centroid(area = True)
            x = x + cent[0] * ar
            y = y + cent[1] * ar
            a = a + ar
        x = x / a
        y = y / a
        if area:
            return ((x, y), a)
        return (x, y)
    
    def area(self):
        vertices = self.get_vertices()
        if vertices is None:
            return None
        if len(vertices) == 1:
            return vertices[0].circles[0].area()
        if len(vertices) == 2:
            s_seg = min(vertices[0].circles).to_segment(vertices[0],vertices[1])
            l_seg = max(vertices[0].circles).to_segment(vertices[0],vertices[1])
            s_a = s_seg.area()
            l_a = l_seg.inv_area()
            return s_a + l_a
        shapes = self.get_circle_segments()
        shapes.append(Polygon(vertices))
        return sum(map(lambda x: x.area(), shapes))
    
    def get_circle_segments(self):
        if self.get_vertices() is None or len(self.get_vertices()) == 1:
            return None
        vertices = Polygon(self.get_vertices()).vertices.copy()
        vertices.append(vertices[0])
        segments = []
        for i in range(0,len(vertices)-1):
            circle = self.matching_circle(vertices[i], vertices[i+1])
            segments.append(circle.to_segment(vertices[i],vertices[i+1]))
        return segments
    
    def matching_circle(self, vertice1, vertice2):
        circles = []
        for c1 in vertice1.circles:
            if c1 in vertice2.circles:
                circles.append(c1)
        return max(circles)
    
    def poly_centroid(self):
        return Polygon(self.vertices).centroid()
    def average(self):
        return Polygon(self.vertices).average()

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
        return [Intersection((x1,y1),(self,other)), \
                Intersection((x2,y2),(self,other))]
        
    def contains(self, pt):
        return self.r > dist(self.coord, pt.coord)
    
    def to_point(self):
        return Intersection(self.coord, [self])
    
    def to_segment(self, pt1, pt2):
        return CircleSegment(self, pt1.angle(self.coord), pt2.angle(self.coord))

    def area(self):
        return pi*self.r**2
    
    def __hash__(self):
        return hash(hash(self.coord) + hash(self.r))
    def __eq__(self, other):
        return self.coord == other.coord and self.r == other.r
    def __gt__(self, other):
        return self.r > other.r

class Point:
    def __init__(self, coord):
        self.coord = coord
    
    def angle(self, centre):
        theta = atan2(self.coord[1]-centre[1], self.coord[0]-centre[0])
        if theta < 0:
            return theta + 2*pi
        return theta
    
    def __hash__(self):
        return hash(self.coord)
    def __eq__(self, other):
        return self.coord == other.coord
    def __str__(self):
        return str(self.coord)
    def __repr__(self):
        return self.__str__()

class Intersection(Point):
    def __init__(self, coord, circles):
        super(Intersection, self).__init__(coord)
        self.circles = circles

class Polygon:
    def __init__(self, vertices):
        self.vertices = vertices
        self.order_vertices()
    
    def centroid(self, area = False):
        if len(self.vertices) < 3:
            return self.average()
        vertices = self.vertices.copy()
        vertices.append(vertices[0])
        a, x, y = (0, 0, 0)
        for i in range(0,len(vertices)-1):
            temp_a = (vertices[i].coord[0]*vertices[i+1].coord[1] - \
                     vertices[i+1].coord[0]*vertices[i].coord[1])
            x = x + (vertices[i].coord[0] + vertices[i+1].coord[0]) * temp_a
            y = y + (vertices[i].coord[1] + vertices[i+1].coord[1]) * temp_a
            a = a + temp_a
        a = a/2
        c = (x/(6*a),y/(6*a))
        if area:
            return (c, a)
        return c
    
    def area(self):
        if len(self.vertices) < 3:
            return 0
        vertices = self.vertices.copy()
        vertices.append(vertices[0])
        a = 0
        for i in range(0,len(vertices)-1):
            a = a + vertices[i].coord[0]*vertices[i+1].coord[1] - \
                    vertices[i+1].coord[0]*vertices[i].coord[1]
        return a / 2
    
    def order_vertices(self):
        centre = self.average()
        self.vertices.sort(key = lambda x: x.angle(centre))
        return self.vertices
    
    def average(self):
        if self.vertices is False:
            return None
        return (mean(list(map(lambda x: x.coord[0], self.vertices))), \
                mean(list(map(lambda x: x.coord[1], self.vertices))))
        
class CircleSegment:
    def __init__(self, circle, psi1, psi2):
        self.coord = circle.coord
        self.r = circle.r
        self.psi1 = psi1
        self.psi2 = psi2
    
    def centroid(self, area = False):
        theta = self.angle()
        d = (4*self.r*sin(theta/2)**3)/(3*(theta-sin(theta)))
        theta_c = self.centre_angle()
        c = (d*cos(theta_c)+self.coord[0], d*sin(theta_c)+self.coord[1])
        if area:
            return (c, self.area())
        return c
        
    def area(self):
        theta = self.angle()
        return (self.r**2/2)*(theta-sin(theta))
        
    def angle(self):       
        theta = abs(self.psi1 - self.psi2)
        return min(theta, 2*pi - theta)
    
    def centre_angle(self):
        if 0 < self.psi2 - self.psi1 < pi or \
           (pi <= self.psi1 and self.psi2 < (self.psi1+pi) % (2*pi)):
            return (self.psi1+self.angle()/2) % (2*pi)
        return (self.psi2+self.angle()/2) % (2*pi)
    
    def inv_centroid(self, area = False):
        theta = 2*pi - self.angle()
        d = (2*self.r*sin(theta/2))/(3*theta/2)
        theta_c = (self.centre_angle()+pi) % (2*pi)
        sector_c = (d*cos(theta_c)+self.coord[0], d*sin(theta_c)+self.coord[1])
        sector_a = theta*self.r**2
        tri_c = (self.coord[0]+(self.r*cos(self.psi1)+self.r*cos(self.psi2))/3, \
                 self.coord[1]+(self.r*sin(self.psi1)+self.r*sin(self.psi2))/3)
        tri_a = self.r**2 * sin(self.angle()) / 2
        a = sector_a + tri_a
        x = (sector_c[0]*sector_a + tri_c[0]*tri_a) / a
        y = (sector_c[1]*sector_a + tri_c[1]*tri_a) / a
        if area:
            return ((x,y), a)
        return (x,y)
        
    def inv_area(self):
        return Circle(self.coord, self.r).area() - self.area()
    
    def inv_contains(self, pt):
        d = self.r*cos(self.angle()/2)
        return d > dist(self.coord, pt.coord)
        
