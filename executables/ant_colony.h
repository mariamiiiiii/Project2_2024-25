#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <stdexcept> // For runtime_error
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>

// Define CGAL kernel and types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

int ant_colony(std::vector<Point> points, DT& dt, int L, int Kappa, int max_iterations);