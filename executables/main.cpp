#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Constrained_Delaunay_triangulation_2.h>
#include <CGAL/draw_triangulation_2.h>
#include <vector>
#include <utility>
#include <iostream>
#include "inputs.h"
#include "flipEdges.h"
#include "output.h"
#include "centroid.h"
#include "projection.h"
#include "center.h"
#include "circumcenter.h"
#include "inside_convex_polygon_centroid.h"

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> CDT;
typedef CDT::Point Point;


using namespace std;

int main(int argc, char* argv[]) {
    // Initialize the Constrained Delaunay Triangulation (CDT)
    CDT cdt;

    // Get data from the executable function
    InputData input = inputs();

    // Get points
    vector<Point> points = input.points;

    // Get Region Boundary
    vector<int> region_boundary = input.region_boundary;

    // Insert points into the triangulation
    for (const Point& p : points) {
        cdt.insert(p);
    }

    // Insert the region boundary as a constrained polygon
    std::vector<Point> polygon;
    for (int idx : region_boundary) {
        if (idx < points.size()) {
            polygon.push_back(points[idx]);
        } else {
            cerr << "Invalid index in region_boundary: " << idx << endl;
        }
    }

    // Append the first point again to close the polygon
    if (!region_boundary.empty()) {
        int first_idx = region_boundary[0];
        if (first_idx < points.size()) {
            polygon.push_back(points[first_idx]);
        } else {
            cerr << "Invalid first index in region_boundary: " << first_idx << endl;
        }
    }

    // Check if the polygon is valid and insert the constraint
    if (polygon.size() > 2) {
        cdt.insert_constraint(polygon.begin(), polygon.end());
    } else {
        cerr << "Not enough points to form a boundary." << endl;
    }

    // Define and add the constrained edges (from additional_constraints)
    const std::vector<std::vector<int>>& constraints = input.additional_constraints;

    // Insert constrained edges based on the provided indices
    for (const auto& constraint : constraints) {  
        if (constraint.size() == 2) {
            int idx1 = constraint[0];
            int idx2 = constraint[1];
            if (idx1 < points.size() && idx2 < points.size()) {
                cdt.insert_constraint(points[idx1], points[idx2]);
            } else {
                cerr << "Invalid constraint index: " << idx1 << ", " << idx2 << endl;
            }
        }
    }

    // Prompt user to choose the Steiner point insertion method
    cout << "Please choose a method for Steiner points from the following options:\n";
    cout << "1: Center of longest edge\n";
    cout << "2: Projection\n";
    cout << "3: Circumcenter\n";
    cout << "4: Centroid of internal convex polygon\n";
    cout << "5: Centroid\n";
    cout << "6: Flip\n";
    cout << "Enter the number corresponding to your choice: ";

    int choice;
    cin >> choice;

    // Execute the chosen method based on user input
    switch (choice) {
        case 1:
            center_steiner_points(points, cdt);
            break;
        case 2:
            projection(points, cdt);
            break;
        case 3:
            circumcenter_steiner_points(points, cdt);
            break;
        case 4:
            inside_convex_polygon_centroid_steiner_points(points, cdt);
            break;
        case 5:
            centroid_steiner_points(points, cdt);
            break;
        case 6:
            flip_edges(points, cdt);
            break;
        default:
            cerr << "Invalid choice. Please enter a number between 1 and 5.\n";
            return 1;
    }

    return 0;
}
