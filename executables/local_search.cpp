#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <cmath>
#include "output.h"
#include "utils.h"
#include "local_search.h"
#include "projection.h" // Assuming projection method is defined here
#include "circumcenter.h" // Assuming circumcenter method is defined here
#include "centroid.h" // Assuming centroid method is defined here

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Point Point;
typedef DT::Edge Edge;
typedef DT::Face_handle FaceHandle;

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