#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <string>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;

void output(const std::vector<std::pair<Point, Point>>& edges, std::vector<Point> steiner_points_given);
