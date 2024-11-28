#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <cmath>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;


double calculateEnergy(DT& dt, double alpha, double beta);
// Function to perform Simulated Annealing
std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner_sa(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points);
int simulated_annealing(std::vector<Point> points, DT dt, double alpha, double beta, int max_iterations);

