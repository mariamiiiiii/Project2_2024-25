#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <string>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;

Point project_point_onto_line(const Point& P, const Point& A, const Point& B);

DT projection(std::vector<Point> points, DT dt);