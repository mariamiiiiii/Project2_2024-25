#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
#include "output.h"
#include "local_search.h"
#include "projection.h" // Assuming projection method is defined here
#include "circumcenter.h" // Assuming circumcenter method is defined here
#include "centroid.h" // Assuming centroid method is defined here

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

// Local search function
std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points) {
    bool added_steiner = false;

        // Store obtuse triangles
        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            int obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                // Get the vertices of the obtuse triangle
                Point p_obtuse = face->vertex(obtuse_vertex)->point();
                Point p1 = face->vertex((obtuse_vertex + 1) % 3)->point();
                Point p2 = face->vertex((obtuse_vertex + 2) % 3)->point();

                // Test multiple methods for Steiner point placement
                Point projection_point = project_point_onto_line(p_obtuse, p1, p2);
                Point circumcenter_point = circumcenter(p_obtuse, p1, p2);
                Point centroid_point = calculate_centroid(p_obtuse, p1, p2);

                // Evaluate obtuse triangle reduction for each method
                DT temp_dt = dt;
                temp_dt.insert(projection_point);
                int count_projection = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(circumcenter_point);
                int count_circumcenter = count_obtuse_triangles(temp_dt);

                temp_dt = dt;
                temp_dt.insert(centroid_point);
                int count_centroid = count_obtuse_triangles(temp_dt);

                // Select the best point
                Point best_point;
                int min_count = std::min({count_projection, count_circumcenter, count_centroid});
                if (min_count == count_projection) {
                    best_point = projection_point;
                } else if (min_count == count_circumcenter) {
                    best_point = circumcenter_point;
                } else {
                    best_point = centroid_point;
                }

                // Add the best point to the DT
                points.push_back(best_point);
                steiner_points.push_back(best_point);
                added_steiner = true;
            }
        }

    for (const Point& p : steiner_points) {
        dt.insert(p);
    }

    return {steiner_points, points};
}

int local_search(std::vector<Point> points, DT dt, int max_iterations) {
    bool obtuse_exists = true;
    int obtuse_count = 0;
    int iterations = 0;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<std::pair<typename DT::Point, typename DT::Point>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    
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

    edges = print_edges(dt);
    output(edges, steiner_points);
    CGAL::draw(dt);

    return 0;
}