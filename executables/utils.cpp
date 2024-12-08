#include <CGAL/Delaunay_triangulation_2.h>  // For Delaunay triangulation
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>  // Kernel
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/Triangle_2.h>              // For working with triangles
#include <cmath>                          // For geometric calculations (e.g., angles)
#include <utility>                        // For std::pair
#include "utils.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

// Function to calculate the angle between two points and a common vertex
template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3) {
    double a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3)));
    double c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

    // Cosine rule to calculate the angle
    double cos_angle = (b * b + c * c - a * a) / (2 * b * c);
    return std::acos(cos_angle) * 180.0 / M_PI;
}

// Function to check if a triangle is obtuse
template <typename FaceHandle>
int obtuse_vertex_index(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    
    if (angle1 > 90.0 + 0.01) return 0;
    if (angle2 > 90.0 + 0.01) return 1;
    if (angle3 > 90.0 + 0.01) return 2;
    return -1;
}

// Helper function to evaluate obtuse angle count
int count_obtuse_triangles(const DT& dt) {
    int count = 0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++count;
        }
    }
    return count;
}

// Function to print the edges of the triangulation
template <typename DT>
std::vector<std::pair<size_t, size_t>> print_edges(const DT& dt, std::vector<Point> points) {
    std::vector<std::pair<size_t, size_t>> edges; // Corrected the type

    // Create a map from Point to its index in the points vector
    std::map<typename DT::Point, size_t> point_index_map;
    for (size_t i = 0; i < points.size(); ++i) {
        point_index_map[points[i]] = i;
    }

    // Print edges and their indices
    for (auto edge = dt.finite_edges_begin(); edge != dt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3)->point();
        auto v2 = edge->first->vertex((edge->second + 2) % 3)->point();
        size_t idx1 = point_index_map[v1];
        size_t idx2 = point_index_map[v2];

        // Store only indices
        edges.emplace_back(idx1, idx2);

    }

    return edges;
}