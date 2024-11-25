#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <vector>
#include <cmath>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;

// Function to perform local search
std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points);
int local_search(std::vector<Point> points, DT dt, int max_iterations);

#endif // LOCAL_SEARCH_H
