#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jul 19 23:25:21 2017

@author: m
"""
from math import sqrt, atan2, sin, cos, pi
from numpy import mean
import sys
import abc
from abc import abstractmethod

if sys.version_info >= (3, 4):
    ABC = abc.ABC
else:
    ABC = abc.ABCMeta('ABC', (), {})


class Shape(ABC):
    @abstractmethod
    def centroid(self, area = False):
        pass    
    @abstractmethod
    def area(self):
        pass  
    @abstractmethod
    def contains_point(self, pt):
        pass

class Point:
    def __init__(self, coord):
        self.x = coord[0]
        self.y = coord[1]
    
    def polar_angle(self, origin):
        theta = atan2(self.y-origin[1], self.x-origin[0])
        if theta < 0:
            return theta + 2*pi
        return theta   
    def coord(self):
        return (self.x, self.y)
    def __hash__(self):
        return hash(self.coord())
    def __eq__(self, other):
        return self.coord() == other.coord()
    def __str__(self):
        return str(self.coord())
    def __repr__(self):
        return self.__str__()

class CircleOverlap(Shape):
    def __init__(self, circles):
        self.circles = circles
        self.vertices = None
        
    def centroid(self, area = False):
        vertices = self.get_vertices()
        if vertices is None:
            return None
        if len(vertices) == 1: # if just a circle
            c = vertices[0].coord()
            if area:
                return (c, vertices[0].circles[0].area())
            return c
        if len(vertices) == 2: # if only two intersecting circles
            l_seg = max(vertices[0].circles).get_segment(vertices[0],vertices[1])
            s_seg = min(vertices[0].circles).get_segment(vertices[0],vertices[1])
            if l_seg.get_major_segment().contains_point(s_seg):
                s_seg = s_seg.get_major_segment()
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
            self.vertices = [circle]
            return
        if len(self.vertices) == 1: # only a circle
            intersection = self.vertices[0].intersect(circle)
            if intersection is None:
                self.vertices = False
                return
            if type(intersection) is Circle:
                self.vertices = [intersection]
                return
            self.vertices = intersection
            return
        circles = self.filter_vertices(circle)
        if len(circles) == 0: # circle within vertices
            circles = self.get_defining_circles()
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
            if circle.contains_point(pt):
                vertices.append(pt)
                circles.add(pt.circles[0])
                circles.add(pt.circles[1])
        if len(vertices) == 0:
            return []
        self.vertices = vertices
        return circles
    
    def get_defining_circles(self):
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
                    if not(c.contains_point(pt)):
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
            vertices = [min(self.circles)]
            self.vertices = vertices
            return vertices
        # check for intersections contained within all circles
        vertices = []
        for pt in candidates:
            switch = True
            for c in self.circles:
                if c not in pt.circles and not(c.contains_point(pt)):
                    switch = False
                    break
            if switch:
                vertices.append(pt)
        if len(vertices) == 0: # a circles is contained all within other circles
            vertices = [min(self.circles)]
            self.vertices = vertices
            return vertices
        self.vertices = vertices
        return vertices
      
    def area(self):
        vertices = self.get_vertices()
        if vertices is None:
            return None
        if len(vertices) == 1:
            return vertices[0].circles[0].area()
        if len(vertices) == 2:
            l_seg = max(vertices[0].circles).get_segment(vertices[0],vertices[1])
            s_seg = min(vertices[0].circles).get_segment(vertices[0],vertices[1])
            if l_seg.get_major_segment().contains_point(s_seg):
                s_seg = s_seg.get_major_segment()
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
            segments.append(circle.get_segment(vertices[i],vertices[i+1]))
        return segments
    
    def matching_circle(self, vertice1, vertice2):
        circles = []
        for c1 in vertice1.circles:
            if c1 in vertice2.circles:
                circles.append(c1)
        return max(circles)
    
    def contains_point(self, pt):
        circles = self.get_defining_circles()
        for c in circles:
            if not(c.contains_point(pt)):
                return False
        return True

class Circle(Shape, Point):
    def __init__(self, coord, r):
        Point.__init__(self, coord)
        self.r = float(r)
    
    def centroid(self, area = False):
        if area:
            return (self.coord(), self.area())
        return self.coord()
    
    def intersect(self, other):
        d = dist(self.coord(), other.coord())
        if d >= self.r + other.r or (d == 0 and self.r == other.r):
            return None
        if d <= abs(self.r - other.r):
            return min(self, other)
        l = (self.r**2 - other.r**2 + d**2) / (2*d)
        h = sqrt(self.r**2 - l**2)
        xmid = self.x + l/d*(other.x - self.x)
        ymid = self.y + l/d*(other.y - self.y)
        x1 = xmid + h/d*(other.y - self.y)
        y1 = ymid - h/d*(other.x - self.x)
        x2 = xmid - h/d*(other.y - self.y)
        y2 = ymid + h/d*(other.x - self.x)
        return [Intersection((x1,y1),(self,other)), \
                Intersection((x2,y2),(self,other))]
        
    def contains_point(self, pt):
        return self.r > dist(self.coord(), pt.coord())
       
    def get_segment(self, pt1, pt2, minor = True):
        psi1 = pt1.polar_angle(self.coord())
        psi2 = pt2.polar_angle(self.coord())
        if minor:
            return CircleSegment(self.coord(), self.r, (psi1, psi2))
        return MajorSegment(self.coord(), self.r, (psi1, psi2))
    
    def area(self):
        return pi*self.r**2
    
    def __hash__(self):
        return hash(hash(self.coord()) + hash(self.r))
    def __eq__(self, other):
        return self.coord == other.coord() and self.r == other.r
    def __gt__(self, other):
        return self.r > other.r

class Intersection(Point):
    def __init__(self, coord, circles):
        super(Intersection, self).__init__(coord)
        self.circles = circles

class Polygon(Shape):
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
            temp_a = (vertices[i].x*vertices[i+1].y - \
                     vertices[i+1].x*vertices[i].y)
            x = x + (vertices[i].x + vertices[i+1].x) * temp_a
            y = y + (vertices[i].y + vertices[i+1].y) * temp_a
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
            a = a + vertices[i].x*vertices[i+1].y - \
                    vertices[i+1].x*vertices[i].y
        return a / 2
    
    def contains_point(self, pt):
        return False # not implemented yet
    
    def order_vertices(self):
        centre = self.average()
        self.vertices.sort(key = lambda x: x.polar_angle(centre))
        return self.vertices
    
    def average(self):
        if self.vertices is False:
            return None
        return (mean(list(map(lambda x: x.x, self.vertices))), \
                mean(list(map(lambda x: x.y, self.vertices))))
        
class CircleSegment(Circle):
    def __init__(self, coord, r, polar_angles):
        super(CircleSegment, self).__init__(coord, r)
        if 0 < polar_angles[1] - polar_angles[0] < pi or \
           (pi <= polar_angles[0] and polar_angles[1] < (polar_angles[0]+pi) % (2*pi)):
               self.psi1 = polar_angles[0]
               self.psi2 = polar_angles[1]
        else:
            self.psi1 = polar_angles[1]
            self.psi2 = polar_angles[0]
    
    def centroid(self, area = False):
        theta = self.angle()
        d = (4*self.r*sin(theta/2)**3)/(3*(theta-sin(theta)))
        theta_c = self.centre_polar_angle()
        c = (d*cos(theta_c)+self.x, d*sin(theta_c)+self.y)
        if area:
            return (c, self.area())
        return c
        
    def area(self):
        theta = self.angle()
        return (self.r**2/2)*(theta-sin(theta))
        
    def angle(self):       
        theta = abs(self.psi2 - self.psi1)
        return min(theta, 2*pi - theta)
    
    def centre_polar_angle(self):
        return (self.psi1+self.angle()/2) % (2*pi)
    
    def contains_point(self, pt):
        theta = pt.polar_angle(self.coord())
        return super(CircleSegment, self).contains_point(pt) and \
            self.within_angle(theta)
    
    def within_angle(self, theta):
        return 0 < theta - self.psi1 < self.angle() or \
            0 < theta + 2*pi - self.psi1 < self.angle()
        
    def get_major_segment(self):
        return MajorSegment(self.coord(), self.r, (self.psi1, self.psi2))
    
class MajorSegment(CircleSegment):
    def __init__(self, coord, r, polar_angles):
        super(MajorSegment, self).__init__(coord, r, polar_angles)
        
    def centroid(self, area = False):
        theta = 2*pi - self.angle()
        d = (2*self.r*sin(theta/2))/(3*theta/2)
        theta_c = (self.centre_polar_angle()+pi) % (2*pi)
        sector_c = (d*cos(theta_c)+self.x, d*sin(theta_c)+self.y)
        sector_a = theta*self.r**2
        tri_c, tri_a = self.get_triangle().centroid(area = True)
        a = sector_a + tri_a
        x = (sector_c[0]*sector_a + tri_c[0]*tri_a) / a
        y = (sector_c[1]*sector_a + tri_c[1]*tri_a) / a
        if area:
            return ((x,y), a)
        return (x,y)
    
    def area(self):
        return super(CircleSegment, self).area() - super(MajorSegment, self).area()
    
    def contains_point(self, pt):
        theta = pt.polar_angle(self.coord())
        if self.within_angle(theta):
            return self.get_triangle().contains_point(pt)
        return super(CircleSegment, self).contains_point(pt)
    
    def get_triangle(self):
        vertex1 = Point((self.r*cos(self.psi1)+self.x, self.r*sin(self.psi1)+self.y))
        vertex2 = Point((self.r*cos(self.psi2)+self.x, self.r*sin(self.psi2)+self.y))
        return Triangle([self, vertex1, vertex2])
        
class Triangle(Polygon):
    def __init__(self, vertices):
        if len(vertices) != 3:
            raise NameError('Incorrect number of vertices')
        super(Triangle, self).__init__(vertices)
             
    def centroid(self, area = False):
        c = self.average()
        if area:
            return (c, self.area())
        return c
    
    def area(self):
        return abs(self.signed_area())
    
    def contains_point(self, pt):
        x, y = pt.coord()
        x1, y1 = self.vertices[0].coord()
        x2, y2 = self.vertices[1].coord()
        x3, y3 = self.vertices[2].coord()
        a = (x2*y3-x3*y2+(y2-y3)*x+(x3-x2)*y)/(2*self.signed_area())
        b = (x3*y1-x1*y3+(y3-y1)*x+(x1-x3)*y)/(2*self.signed_area())
        return (0 <= a and 0 <= b and (a + b <= 1))

    def signed_area(self):
        x1, y1 = self.vertices[0].coord()
        x2, y2 = self.vertices[1].coord()
        x3, y3 = self.vertices[2].coord()
        return (x1*y2+x2*y3+x3*y1-x1*y3-x3*y2-x2*y1)/2

def dist(coord1, coord2):
    return sqrt((coord1[0]-coord2[0])**2 + (coord1[1]-coord2[1])**2)
