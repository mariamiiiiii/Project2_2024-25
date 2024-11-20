#ifndef LOCAL_SEARCH_H
#define LOCAL_SEARCH_H

#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <vector>
#include <cmath>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;

// Function to perform local search
std::vector<Point> local_search(CDT& cdt, int max_iterations);

#endif // LOCAL_SEARCH_H
