#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Polygon_2.h>
#include <CGAL/convex_hull_2.h>
#include <vector>
#include <set>
#include <stack>
#include <cmath> // For angle calculations
#include "inside_convex_polygon_centroid.h"
#include "output.h"

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;  // Constrained Delaunay triangulation
typedef K::Point_2 Point;
typedef CGAL::Polygon_2<K> Polygon_2;
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
    return std::acos(cos_angle) * 180.0 / M_PI; // Return angle in degrees
}

// Function to check if a triangle is obtuse
// Returns the index of the obtuse angle's vertex or -1 if no obtuse angle is found
template <typename FaceHandle>
int obtuse_vertex_index(const FaceHandle& face) {
    double angle1 = angle_between(face->vertex(0)->point(), face->vertex(1)->point(), face->vertex(2)->point());
    double angle2 = angle_between(face->vertex(1)->point(), face->vertex(2)->point(), face->vertex(0)->point());
    double angle3 = angle_between(face->vertex(2)->point(), face->vertex(0)->point(), face->vertex(1)->point());
    
    if (angle1 > 90.0 + 0.01) return 0;
    if (angle2 > 90.0 + 0.01) return 1;
    if (angle3 > 90.0 + 0.01) return 2;

    return -1; // No obtuse angle
}


// Function to compute centroid of a convex polygon
Point compute_centroid(const std::vector<Point>& points) {
    double sum_x = 0.0, sum_y = 0.0;
    for (const auto& point : points) {
        sum_x += CGAL::to_double(point.x());
        sum_y += CGAL::to_double(point.y());
    }
    return Point(sum_x / points.size(), sum_y / points.size());
}

std::vector<Point> find_convex_polygon(DT& dt, FaceHandle start_face) {
    std::set<Point> unique_points; // Use a set to avoid duplicate points
    std::set<FaceHandle> visited_faces; // To keep track of visited faces
    std::stack<FaceHandle> face_stack; // Stack for DFS
    face_stack.push(start_face); // Start from the initial face

    while (!face_stack.empty()) {
        FaceHandle current_face = face_stack.top();
        face_stack.pop();

        // If this face has already been visited, skip it
        if (visited_faces.find(current_face) != visited_faces.end()) {
            continue;
        }

        visited_faces.insert(current_face);

        // Check if the current face is obtuse
        if (obtuse_vertex_index(current_face) != -1) {
            // Add the vertices of the obtuse triangle to the unique points set
            for (int i = 0; i < 3; ++i) {
                Point pt = current_face->vertex(i)->point();
                unique_points.insert(pt);
            }

            // Explore neighboring faces
            for (int i = 0; i < 3; ++i) {
                FaceHandle neighbor_face = current_face->neighbor(i);
                if (!dt.is_infinite(neighbor_face) && visited_faces.find(neighbor_face) == visited_faces.end() && obtuse_vertex_index(neighbor_face) != -1) {
                    face_stack.push(neighbor_face);
                }
            }
        }
    }

    // Create a Polygon_2 with the unique points
    Polygon_2 polygon;
    for (const auto& pt : unique_points) {
        polygon.push_back(pt);
    }

    // Check if the polygon is convex
    if (polygon.is_convex()) {
        return std::vector<Point>(polygon.vertices_begin(), polygon.vertices_end());
    } else {
        // If not convex, compute the convex hull
        std::vector<Point> convex_hull;
        CGAL::convex_hull_2(unique_points.begin(), unique_points.end(), std::back_inserter(convex_hull));
        return convex_hull;
    }
}

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

// Function to add Steiner points at the center of convex polygons of obtuse triangles
template <typename DT>
std::vector<Point> add_steiner_in_convex_polygon_centroid(DT& dt, std::vector<Point> steiner_points) {
    bool added_steiner = false;
    int count = 0;

    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        count++;
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            auto polygon_points = find_convex_polygon(dt, face);
            if (!polygon_points.empty()) {  // Only proceed if convex polygon found
                Point centroid = compute_centroid(polygon_points);
                steiner_points.push_back(centroid);
                added_steiner = true;
            }
        }
    }

    for (const Point& p : steiner_points) {
        dt.insert(p);
    }

    return steiner_points;
}

int inside_convex_polygon_centroid_steiner_points(std::vector<Point> points, DT dt) {
    std::vector<Point> steiner_points;
    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;

    bool obtuse_exists = true;
    int obtuse_count = 0;
    int iterations = 0;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    while (obtuse_exists && iterations <= 5) {
        
        steiner_points = add_steiner_in_convex_polygon_centroid(dt, steiner_points);
        
        obtuse_exists = false;
        obtuse_count = 0;  

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            int obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
                obtuse_count++;
                CGAL::draw(dt);
            }
        }
        iterations++;
    }

    edges = print_edges(dt);

    output(edges, steiner_points);
    CGAL::draw(dt);

    return 0;
    
}
