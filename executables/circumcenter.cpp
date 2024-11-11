#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <CGAL/Polygon_2.h>
#include <CGAL/Polygon_2_algorithms.h>
#include "output.h"

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT; 
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;
typedef K::FT FT;
typedef CGAL::Polygon_2<K> Polygon;


// Συνάρτηση για έλεγχο αν ένα σημείο είναι εντός του κυρτού περιβλήματος
bool is_within_convex_hull(const Point& p, const std::vector<Point>& convex_hull) {
    Polygon hull_polygon(convex_hull.begin(), convex_hull.end());
    return CGAL::bounded_side_2(hull_polygon.vertices_begin(), hull_polygon.vertices_end(), p, K()) 
           != CGAL::ON_UNBOUNDED_SIDE;
}


// Function to calculate the angle between two points and a common vertex
template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3) {
    double a = std::sqrt(CGAL::to_double(CGAL::squared_distance(p2, p3)));
    double b = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p3)));
    double c = std::sqrt(CGAL::to_double(CGAL::squared_distance(p1, p2)));

    // Cosine rule to calculate the angle
    double cos_angle = (b * b + c * c - a * a) / (2 * b * c);
    return std::acos(cos_angle) * 180.0 / M_PI; // Return angle in degrees
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

template <typename FaceHandle>
std::pair<int, double> obtuse_vertex_index_and_angle(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    
    if (angle1 > 90.0) return std::make_pair(0, angle1);
    if (angle2 > 90.0) return std::make_pair(1, angle2);
    if (angle3 > 90.0) return std::make_pair(2, angle3);
    
    return std::make_pair(-1, 0.0);  // No obtuse angle
}

// Function to calculate the circumcenter of a triangle
Point circumcenter(const Point& p1, const Point& p2, const Point& p3) {
    // Using the CGAL function to calculate the circumcenter
    return CGAL::circumcenter(p1, p2, p3);
}



// Function to print the edges of the triangulation
template <typename DT>
std::vector<std::pair<typename DT::Point, typename DT::Point>> print_edges(const DT& dt) {
    // Define a vector to hold pairs of points representing edges
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;
    for (auto edge = dt.finite_edges_begin(); edge != dt.finite_edges_end(); ++edge) {
        auto v1 = edge->first->vertex((edge->second + 1) % 3)->point();
        auto v2 = edge->first->vertex((edge->second + 2) % 3)->point();
        // Add the edge to the vector
        edges.emplace_back(v1, v2);
    }

    // Return the vector of edges
    return edges;
}

template <typename DT>
std::vector<Point> add_steiner_in_circumcenter(DT& dt, std::vector<Point> steiner_points,  const std::vector<Point>& convex_hull) {
    bool added_steiner = false;

    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            Point p1 = face->vertex(0)->point();
            Point p2 = face->vertex(1)->point();
            Point p3 = face->vertex(2)->point();
            Point circumcenter_point = circumcenter(p1, p2, p3);

            // Έλεγχος αν το Steiner σημείο βρίσκεται εντός του κυρτού περιβλήματος
            if (is_within_convex_hull(circumcenter_point, convex_hull)) {
                steiner_points.push_back(circumcenter_point);
                added_steiner = true;
            }
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

int circumcenter_steiner_points(std::vector<Point> points, DT dt) {

    std::vector<Point> convex_hull;
    CGAL::convex_hull_2(points.begin(), points.end(), std::back_inserter(convex_hull));
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;

    bool obtuse_exists = true;
    int obtuse_count = 0;
    int iterations = 0;
    
    std::vector<Point> steiner_points;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    obtuse_count = 0;  
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            obtuse_count++;
        }
    }

    while (obtuse_exists && iterations <= 5) {
        
        steiner_points = add_steiner_in_circumcenter(dt, steiner_points, convex_hull);
        
        obtuse_exists = false;
        obtuse_count = 0;

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            auto obtuse = obtuse_vertex_index_and_angle(face);
            auto obtuse_vertex = std::get<0>(obtuse);
            auto obtuse_angle = std::get<1>(obtuse);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
                obtuse_count++;

                // Get coordinates of the obtuse vertex and print detailed information
                Point obtuse_point = face->vertex(obtuse_vertex)->point();

                // Print the coordinates of the other two vertices
                Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
                Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();
            }
        }
        iterations++;
    }

    edges = print_edges(dt);

    output(edges, steiner_points);
    CGAL::draw(dt);

    return 0;
}
