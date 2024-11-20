#include "local_search.h"
#include "projection.h" // Assuming projection method is defined here
#include "circumcenter.h" // Assuming circumcenter method is defined here
#include "centroid.h" // Assuming centroid method is defined here

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
int count_obtuse_triangles(const CDT& cdt) {
    int count = 0;
    for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++count;
        }
    }
    return count;
}

// Local search function
std::vector<Point> local_search(CDT& cdt, int max_iterations) {
    std::vector<Point> all_steiner_points;
    bool improved = true;
    int iteration = 0;

    while (improved && iteration < max_iterations) {
        improved = false;

        // Store obtuse triangles
        for (auto face = cdt.finite_faces_begin(); face != cdt.finite_faces_end(); ++face) {
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
                CDT temp_cdt = cdt;
                temp_cdt.insert(projection_point);
                int count_projection = count_obtuse_triangles(temp_cdt);

                temp_cdt = cdt;
                temp_cdt.insert(circumcenter_point);
                int count_circumcenter = count_obtuse_triangles(temp_cdt);

                temp_cdt = cdt;
                temp_cdt.insert(centroid_point);
                int count_centroid = count_obtuse_triangles(temp_cdt);

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

                // Add the best point to the CDT
                cdt.insert(best_point);
                all_steiner_points.push_back(best_point);
                improved = true;
            }
        }

        ++iteration;
    }

    return all_steiner_points;
}
