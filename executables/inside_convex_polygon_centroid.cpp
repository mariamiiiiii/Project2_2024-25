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

std::vector<Point> find_convex_polygon(DT& dt, FaceHandle current_face) {
    std::vector<Point> steiner_points;
    std::vector<Point> polygon_points;

    // Helper lambda to check convexity and add a Steiner point
    auto process_polygon = [&](const std::vector<Point>& points) -> bool {
        Polygon_2 polygon(points.begin(), points.end());
        if (!polygon.is_convex()) {
            return false;
        }
        Point centroid = compute_centroid(points);
        steiner_points.push_back(centroid);
        return true;
    };

    // Collect points of the current face
    for (int i = 0; i < 3; ++i) {
        if (!dt.is_infinite(current_face->vertex(i))) {
            polygon_points.push_back(current_face->vertex(i)->point());
        }
    }

    // Check if the current face and all its neighbors form a convex polygon
    std::vector<FaceHandle> neighbors;
    for (int i = 0; i < 3; ++i) {
        FaceHandle neighbor = current_face->neighbor(i);
        if (!dt.is_infinite(neighbor)) {
            neighbors.push_back(neighbor);
        }
    }

    if (neighbors.size() == 3) {
        for (const auto& neighbor : neighbors) {
            for (int i = 0; i < 3; ++i) {
                if (!dt.is_infinite(neighbor->vertex(i))) {
                    Point p = neighbor->vertex(i)->point();
                    if (std::find(polygon_points.begin(), polygon_points.end(), p) == polygon_points.end()) {
                        polygon_points.push_back(p);
                    }
                }
            }
        }

        if (process_polygon(polygon_points)) {
            return steiner_points;
        }
    }

    // Try current face with pairs of neighbors
    for (size_t i = 0; i < neighbors.size(); ++i) {
        for (size_t j = i + 1; j < neighbors.size(); ++j) {
            polygon_points.clear();

            // Add current face points
            for (int k = 0; k < 3; ++k) {
                if (!dt.is_infinite(current_face->vertex(k))) {
                    polygon_points.push_back(current_face->vertex(k)->point());
                }
            }

            // Add points from the two neighbors
            for (const auto& neighbor : {neighbors[i], neighbors[j]}) {
                for (int k = 0; k < 3; ++k) {
                    if (!dt.is_infinite(neighbor->vertex(k))) {
                        Point p = neighbor->vertex(k)->point();
                        if (std::find(polygon_points.begin(), polygon_points.end(), p) == polygon_points.end()) {
                            polygon_points.push_back(p);
                        }
                    }
                }
            }

            if (process_polygon(polygon_points)) {
                return steiner_points;
            }
        }
    }

    // Fallback: Current face with one neighbor
    for (const auto& neighbor : neighbors) {
        polygon_points.clear();

        // Add current face points
        for (int k = 0; k < 3; ++k) {
            if (!dt.is_infinite(current_face->vertex(k))) {
                polygon_points.push_back(current_face->vertex(k)->point());
            }
        }

        // Add points from the neighbor
        for (int k = 0; k < 3; ++k) {
            if (!dt.is_infinite(neighbor->vertex(k))) {
                Point p = neighbor->vertex(k)->point();
                if (std::find(polygon_points.begin(), polygon_points.end(), p) == polygon_points.end()) {
                    polygon_points.push_back(p);
                }
            }
        }

        if (process_polygon(polygon_points)) {
            return steiner_points;
        }
    }

    return steiner_points;
}


// std::vector<Point> find_convex_polygon(DT& dt, FaceHandle start_face) {
//     std::set<Point> unique_points; // Use a set to avoid duplicate points
//     std::set<FaceHandle> visited_faces; // To keep track of visited faces
//     std::stack<FaceHandle> face_stack; // Stack for DFS
//     face_stack.push(start_face); // Start from the initial face

//     while (!face_stack.empty()) {
//         FaceHandle current_face = face_stack.top();
//         face_stack.pop();

//         // If this face has already been visited, skip it
//         if (visited_faces.find(current_face) != visited_faces.end()) {
//             continue;
//         }

//         visited_faces.insert(current_face);

//         // Check if the current face is obtuse
//         if (obtuse_vertex_index(current_face) != -1) {
//             // Add the vertices of the obtuse triangle to the unique points set
//             for (int i = 0; i < 3; ++i) {
//                 if (!dt.is_infinite(current_face->vertex(i))) {
//                     unique_points.insert(current_face->vertex(i)->point());
//                 }
//             }

//             // Explore neighboring faces
//             for (int i = 0; i < 3; ++i) {
//                 FaceHandle neighbor_face = current_face->neighbor(i);
//                 if (!dt.is_infinite(neighbor_face) && visited_faces.find(neighbor_face) == visited_faces.end() && obtuse_vertex_index(neighbor_face) != -1) {
//                     face_stack.push(neighbor_face);
//                 }
//             }
//         }
//     }

//     // Create a Polygon_2 with the unique points
//     Polygon_2 polygon;
//     for (const auto& pt : unique_points) {
//         polygon.push_back(pt);
//     }

//     // Check if the polygon is convex
//     if (polygon.is_convex()) {
//         return std::vector<Point>(polygon.vertices_begin(), polygon.vertices_end());
//     } else {
//         // If not convex, compute the convex hull
//         std::vector<Point> convex_hull;
//         CGAL::convex_hull_2(unique_points.begin(), unique_points.end(), std::back_inserter(convex_hull));
//         return convex_hull;
//     }
// }

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

        std::cout << "Edge between indices " << idx1 << " and " << idx2 << std::endl;
    }

    // Print point indices for all vertices in the triangulation
    for (auto vertex = dt.finite_vertices_begin(); vertex != dt.finite_vertices_end(); ++vertex) {
        auto pt = vertex->point();
        size_t idx = point_index_map[pt];
        std::cout << "Point: " << pt << ", Index: " << idx << std::endl;
    }

    return edges;
}

// Function to add Steiner points at the center of convex polygons of obtuse triangles
template <typename DT>
std::pair<std::vector<Point>, std::vector<Point>> add_steiner_in_convex_polygon_centroid(DT& dt, std::vector<Point> steiner_points, std::vector<Point> points) {

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
                points.push_back(centroid);
                added_steiner = true;
            }
        }
    }

    for (const Point& p : steiner_points) {
        dt.insert(p);
    }

    return {steiner_points, points};
}

int inside_convex_polygon_centroid_steiner_points(std::vector<Point> points, DT dt) {
    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;
    std::vector<std::pair<size_t, size_t>> edges;

    bool obtuse_exists = true;
    int obtuse_count = 0;
    int iterations = 0;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt);

    while (obtuse_exists && iterations <= 5) {
        
        all_points = add_steiner_in_convex_polygon_centroid(dt, steiner_points, points);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        
        obtuse_exists = false;
        obtuse_count = 0;  

        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
            int obtuse_vertex = obtuse_vertex_index(face);
            if (obtuse_vertex != -1) {
                obtuse_exists = true;
                obtuse_count++;
            }
        }
        iterations++;
    }

    edges = print_edges(dt, all_points.first);

    output(edges, steiner_points);
    CGAL::draw(dt);

    return 0;
    
}
