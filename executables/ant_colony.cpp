#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/draw_triangulation_2.h>
// #include <CGAL/sqrt.h> // Ensure the sqrt function is correctly used
#include <vector>
#include <cmath> // For standard math functions like std::sqrt

// Define CGAL kernel and types
typedef CGAL::Exact_predicates_inexact_constructions_kernel K; // Use this kernel for compatibility
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

// Function to compute the height of point P above line segment AB
K::FT height(const Point& P, const Point& A, const Point& B) {
    // Calculate the determinant for the area of the parallelogram
    K::FT det = (B.x() - A.x()) * (A.y() - P.y()) - (A.x() - P.x()) * (B.y() - A.y());
    
    // Length of the line segment AB
    K::FT length_squared = (B.x() - A.x()) * (B.x() - A.x()) + (B.y() - A.y()) * (B.y() - A.y());

    // Compute the height as |det| / sqrt(length_squared)
    if (length_squared == 0) {
        throw std::runtime_error("Segment length is zero; cannot compute height.");
    }
    return CGAL::abs(det) / CGAL::sqrt(length_squared);
}

// Function to calculate the radius-to-height ratio
K::FT radius_to_height_ratio(const FaceHandle& face) {
    // Compute circumradius
    Point A = face->vertex(0)->point();
    Point B = face->vertex(1)->point();
    Point C = face->vertex(2)->point();

    // Calculate circumradius using squared distance
    K::FT circumradius_squared = CGAL::squared_distance(A, B) * CGAL::squared_distance(A, C) * CGAL::squared_distance(B, C)
                                  / (4 * CGAL::abs(CGAL::area(A, B, C)));
    K::FT circumradius = CGAL::sqrt(circumradius_squared);

    // Compute height of triangle base AB
    K::FT h = height(C, A, B);
    
    // return circumradius / h;

    K::FT r = circumradius / h;
    if(r > 2){
        if((r-1)/r>0){
            return (r-1)/r;
        } else {
            return 0;
        }
    } else if (r>=1 && r<=2){
        return r/(2+r);
    } else if (r<1){
        return 1;
    } else return 0;
}

// Dummy implementations for missing functions
std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner(
    DT& dt,
    const std::vector<Point>& steiner_points,
    const std::vector<Point>& points) {
    // Stub: implement Steiner point optimization
    return {steiner_points, points};
}

int obtuse_vertex_index(const FaceHandle& face) {
    // Stub: implement check for obtuse triangle vertices
    return -1;
}

std::vector<std::pair<size_t, size_t>> print_edges(const DT& dt, const std::vector<Point>& points) {
    // Stub: extract edges from the triangulation
    return {};
}

void output(const std::vector<std::pair<size_t, size_t>>& edges, const std::vector<Point>& steiner_points) {
    // Stub: output edges and Steiner points to a file or console
}

// Ant colony optimization function
int ant_colony(std::vector<Point> points, DT& dt, int L, int Kappa, int max_iterations) {
    bool obtuse_exists = true;
    int iterations = 0;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt); // Draw initial triangulation

    for (int cycle = 1; cycle < L; cycle++) {
        for (int ant = 1; ant < Kappa; ant++) {
            // Stub: implement ant behavior
        }
    }
    
    while (obtuse_exists && iterations <= max_iterations) {
        all_points = add_best_steiner(dt, steiner_points, points);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        obtuse_exists = false;

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            if (obtuse_vertex_index(face) != -1) {
                obtuse_exists = true;
                break;
            }
        }
        iterations++;
    }

    edges = print_edges(dt, all_points.first);
    output(edges, steiner_points);
    CGAL::draw(dt); // Draw final triangulation

    return 0;
}
