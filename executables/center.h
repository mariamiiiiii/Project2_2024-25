#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
// typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef Custom_Constrained_Delaunay_triangulation_2<K> CustomDT; // Renamed typedef
typedef CustomDT::Point Point;

Point longest_edge_center(const Point& p1, const Point& p2);

int center_steiner_points(std::vector<Point> points, CustomDT dt);