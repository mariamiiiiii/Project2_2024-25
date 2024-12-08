#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/squared_distance_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <stdexcept>
#include <iostream>
#include <random>
#include <ctime>
#include <cmath>
#include "ant_colony.h"
#include "utils.h"
#include "output.h"

// Define CGAL kernel and types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

// Function to compute the height of point P above line segment AB
K::FT height(const Point& P, const Point& A, const Point& B) {
    // Calculate the determinant for the area of the parallelogram
    K::FT det = (B.x() - A.x()) * (A.y() - P.y()) - (A.x() - P.x()) * (B.y() - A.y());
    
    // Length of the line segment AB
    K::FT length_squared = (B.x() - A.x()) * (B.x() - A.x()) + (B.y() - A.y()) * (B.y() - A.y());

    if (length_squared == K::FT(0)) {
        throw std::runtime_error("Segment length is zero; cannot compute height.");
    }

    double length = CGAL::to_double(length_squared);
    double abs_det = CGAL::to_double(CGAL::abs(det));

    // Return the height as a K::FT
    return K::FT(abs_det / std::sqrt(length));
}

// Function to calculate the radius-to-height ratio
K::FT radius_to_height_ratio(const FaceHandle& face, DT dt) {
    Point A = face->vertex(0)->point();
    Point B = face->vertex(1)->point();
    Point C = face->vertex(2)->point();

    int obtuse_adjacent_triangles_counter = 0;

    if (obtuse_vertex_index(face) != -1) {
            obtuse_adjacent_triangles_counter++;
            // Explore neighboring faces
            for (int i = 0; i < 3; ++i) {
                FaceHandle neighbor_face = face->neighbor(i);
                if (!dt.is_infinite(neighbor_face) && obtuse_vertex_index(neighbor_face) != -1) {
                    obtuse_adjacent_triangles_counter++;
                }
            }
        }

    K::FT area = CGAL::abs(CGAL::area(A, B, C));
    if (area == K::FT(0)) {
        throw std::runtime_error("Degenerate triangle; cannot compute circumradius.");
    }

    K::FT circumradius_squared = CGAL::squared_distance(A, B) * CGAL::squared_distance(A, C) * CGAL::squared_distance(B, C)
                                  / (4 * area);

    double circumradius_double = std::sqrt(CGAL::to_double(circumradius_squared));

    K::FT circumradius = K::FT(circumradius_double);

    K::FT h = height(C, A, B);
    if (h == K::FT(0)) {
        throw std::runtime_error("Height is zero; invalid triangle.");
    }

    K::FT r = circumradius / h;

    //heuristics based on steiner point techniques
    if (obtuse_adjacent_triangles_counter >=2) {
        return K::FT(1);
    } else if (r >= 1 && r <= 2) {
        return r / (2.0 + r);
    } else if (r < 1) {
        return (3.0 - 2.0*r) / 3.0 > K::FT(0) ? (3.0 - 2.0*r) / 3.0 : K::FT(0);
    } else {
        return (r - 1.0) / r > K::FT(0) ? (r - 1.0) / r : K::FT(0);
    }
}

double calculateDelta(DT& dt, double alpha, double beta, int steiner_points_count) {
    int obtuse_count = 0;
    for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
        int obtuse_vertex = obtuse_vertex_index(face);
        if (obtuse_vertex != -1) {
            ++obtuse_count;
        }
    }

    return 1.0/(1 + alpha * obtuse_count + beta * steiner_points_count);
}

double steiner_point_probability(double t, double ) {
    // Initialize random number generator
    static std::mt19937 rng(static_cast<unsigned int>(std::time(nullptr)));
    std::uniform_real_distribution<double> dist(0.0, 1.0); // Uniform distribution [0, 1]

        double probability = std::pow(t, x) * std::pow(h, y);
        // Generate a random number in [0, 1]
        double random = dist(rng);
        // Accept with probability
        return random < probability;

}

std::pair<std::vector<Point>, std::vector<Point>> add_best_steiner(DT& dt, const std::vector<Point>& steiner_points, const std::vector<Point>& points) {
    // Stub: implement Steiner point optimization
    return {steiner_points, points};
}

// Ant colony optimization function
int ant_colony(std::vector<Point> points, DT& dt, int L, int Kappa, int max_iterations, double alpha, double beta, double lamda ) {
    bool obtuse_exists = true;
    int iterations = 0;
    double pheromone_projection = 1.0; 
    double pheromone_circumcenter = 1.0;
    double pheromone_midpoint = 1.0;
    double pheromone_adjacent_obtuse_triangles = 1.0;

    double delta_projection = 0.0;
    double delta_circumcenter = 0.0;
    double delta_midpoint = 0.0;
    double delta_adjacent_obtuse_triangles = 0.0;

    int method_used = 0;

    double sum = 0.0;

    std::vector<Point> steiner_points;
    std::pair<std::vector<Point>, std::vector<Point>> all_points;

    std::vector<double> probabilities;

    std::vector<std::pair<size_t, size_t>> edges;

    // Insert points into the triangulation
    for (const Point& p : points) {
        dt.insert(p);
    }

    CGAL::draw(dt); // Draw initial triangulation

    
    while (obtuse_exists && iterations <= max_iterations) {
        all_points = add_best_steiner(dt, steiner_points, points);
        steiner_points = all_points.first;  // Extract Steiner points
        points = all_points.second;        // Extract updated points
        obtuse_exists = false;

        int obtuse_count = 0, obtuse_previous_count = 0;  //possible na xreiazetai ena gia kathe methodo???

        for (int cycle = 1; cycle < L; cycle++) {
            for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
                auto obtuse_vertex = obtuse_vertex_index(face);
                if (obtuse_vertex != -1) {
                    obtuse_previous_count++;
                }
            }

            for (int ant = 1; ant < Kappa; ant++) {
                for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
                    if (obtuse_vertex_index(face) != -1) {

                        probabilities = calculateProbabilities(tau, eta, chi, psi);

                        //ylopoihsh improve tringulation

                        for (auto face = dt.finite_faces_begin(); face != dt.finite_faces_end(); ++face) {
                            auto obtuse_vertex = obtuse_vertex_index(face);
                            if (obtuse_vertex != -1) {
                                obtuse_count++;
                            }
                        }
                        break;

                        //mexri edo
                    }
                }   
            }
            
            if(method_used == 1){
                if(obtuse_count < obtuse_previous_count){
                    delta_projection = calculateDelta(dt, alpha, beta, steiner_points.size());
                } else {
                    delta_projection = 0.0;
                }
                pheromone_projection = (1 - lamda) * pheromone_projection + delta_projection;
            } else if (method_used == 2) {
                if(obtuse_count < obtuse_previous_count){
                    delta_circumcenter = calculateDelta(dt, alpha, beta, steiner_points.size());
                } else {
                    delta_circumcenter = 0.0;
                }
                pheromone_circumcenter = (1 - lamda) * pheromone_circumcenter + delta_circumcenter;
            } else if (method_used == 3) {
                if(obtuse_count < obtuse_previous_count){
                    delta_midpoint = calculateDelta(dt, alpha, beta, steiner_points.size());
                } else {
                    delta_midpoint = 0.0;
                }
                pheromone_midpoint = (1 - lamda) * pheromone_midpoint + delta_midpoint;
            } else {
                if(obtuse_count < obtuse_previous_count){
                    delta_adjacent_obtuse_triangles = calculateDelta(dt, alpha, beta, steiner_points.size());
                } else {
                    delta_adjacent_obtuse_triangles = 0.0;
                }
                pheromone_adjacent_obtuse_triangles = (1 - lamda) * pheromone_adjacent_obtuse_triangles + delta_adjacent_obtuse_triangles;
            }
        }

        iterations++;
    }

    edges = print_edges(dt, all_points.first);
    output(edges, steiner_points);
    CGAL::draw(dt); // Draw final triangulation

    return 0;
}
