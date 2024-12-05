#include <CGAL/Simple_cartesian.h>
#include <CGAL/squared_distance_2.h>
#include <cmath> // For sqrt function
#include "circumcenter.h"
#include "output.h"

typedef CGAL::Simple_cartesian<double> K; // Change to exact kernel if needed
typedef K::Point_2 Point;
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;


K::FT height(const Point& P, const Point& A, const Point& B) {
    // Calculate the determinant for the area of the parallelogram
    K::FT det = (B.x() - A.x()) * (A.y() - P.y()) - (A.x() - P.x()) * (B.y() - A.y());
    
    // Length of the line segment AB
    K::FT length_squared = (B.x() - A.x()) * (B.x() - A.x()) + (B.y() - A.y()) * (B.y() - A.y());
    
    // Compute the height as |det| / sqrt(length_squared)
    return CGAL::abs(det) / CGAL::sqrt(length_squared);
}

// Function to calculate the radius-to-height ratio
double radius_to_height_ratio(const FaceHandle& face) {
    double circumradius = CGAL::circumradius(face); // Radius of the circumcircle
    double height = height(triangle);
    return circumradius / height;
}


int ant_colony(std::vector<Point> points, DT dt, int L, int K) {
    bool obtuse_exists = true;
    int obtuse_count = 0;
    int iterations = 0;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);


    for (int cycle = 1; cycle < L; cycle++) {
        for (int ant = 1; ant < K; ant++) {

        }
    }
    
    while (obtuse_exists && iterations <= max_iterations) {
        all_points = add_best_steiner(dt, steiner_points, points);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        obtuse_exists = false;
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
            }
        }
        iterations++;
    }

    edges = print_edges(dt, all_points.first);
    output(edges, steiner_points);
    CGAL::draw(dt);

    return 0;
}