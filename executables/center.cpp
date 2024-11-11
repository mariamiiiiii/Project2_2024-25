#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
#include "center.h"
#include "output.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;

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

// Function to check if a triangle is obtuse and return the index and angle of the obtuse vertex
template <typename FaceHandle>
std::pair<int, double> obtuse_vertex_index_and_angle(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    if (angle1 > 90.0 + 0.01) return std::make_pair(0, angle1);
    if (angle2 > 90.0 + 0.01) return std::make_pair(1, angle2);
    if (angle3 > 90.0 + 0.01) return std::make_pair(2, angle3);
    return std::make_pair(-1, 0.0);
}

// Function to get the edges of the triangulation
template <typename DT>
std::vector<std::pair<typename DT::Point, typename DT::Point>> print_edges(const DT& dt) {
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;
    for (auto edge = dt.finite_edges_begin(); edge != dt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3)->point();
        auto v2 = edge->first->vertex((edge->second + 2) % 3)->point();
        edges.emplace_back(v1, v2);
    }
    return edges;
}

// Function to add Steiner points instead of flipping edges when a triangle is obtuse
template <typename DT>
std::vector<Point> add_steiner_if_obtuse_center(DT& dt, std::vector<Point> steiner_points) {
    bool added_steiner = false;
    // std::vector<Point> steiner_points;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            // Get the vertices of the obtuse triangle
            Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
            Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();
            // Calculate the midpoint of the edge opposite the obtuse angle
            K::FT mid_x = (p1.x() + p2.x()) / 2;
            K::FT mid_y = (p1.y() + p2.y()) / 2;
            Point midpoint(mid_x, mid_y);
            // Add the Steiner point to the list
            steiner_points.push_back(midpoint);
            added_steiner = true;
        }
    }
    // Insert Steiner points into the triangulation and re-triangulate
    for (const Point& p : steiner_points) {
        dt.insert(p);
    }
    if (added_steiner) {
        return steiner_points;
    } else {
        return {};
    }
}

int center_steiner_points(std::vector<Point> points, DT dt) {
    std::vector<Point> steiner_points;
    bool obtuse_exists = true;
    int iterations = 0;
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;
    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }
    CGAL::draw(dt);
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
    }
    while (obtuse_exists && iterations <= 5) {
        steiner_points = add_steiner_if_obtuse_center(dt, steiner_points);
        obtuse_exists = false;
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse = obtuse_vertex_index_and_angle(face);
            auto obtuse_vertex = std::get<0>(obtuse);
            auto obtuse_angle = std::get<1>(obtuse);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
            }
        }
        iterations++;
    }
    edges = print_edges(dt);
    output(edges, steiner_points);
    CGAL::draw(dt);
    return 0;
}
