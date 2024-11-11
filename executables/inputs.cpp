#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "inputs.h"

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;

InputData inputs() {
    // Creation of property tree
    boost::property_tree::ptree pt;

    // Read from JSON file
    try {
        read_json("../input.json", pt);
    } catch (const boost::property_tree::json_parser_error &e) {
        std::cerr << "Error reading JSON: " << e.what() << std::endl;
        return {};  // Return empty struct
    }

    // Get data from property tree
    std::string instance_uid = pt.get<std::string>("instance_uid");
    int num_points = pt.get<int>("num_points");

    // Get points
    std::vector<int> points_x;
    for (const auto& point : pt.get_child("points_x")) {
        points_x.push_back(point.second.get_value<int>());
    }

    std::vector<int> points_y;
    for (const auto& point : pt.get_child("points_y")) {
        points_y.push_back(point.second.get_value<int>());
    }

    // Create CGAL Point objects from the x and y coordinates
    std::vector<Point> points;
    for (int i = 0; i < num_points; ++i) {
        points.push_back(Point(points_x[i], points_y[i])); 
    }

    // Retrieve additional constraints
    std::vector<std::vector<int>> additional_constraints;
    for (const auto& row : pt.get_child("additional_constraints")) {
        std::vector<int> row_values;
        for (const auto& value : row.second) {
            row_values.push_back(value.second.get_value<int>());
        }
        additional_constraints.push_back(row_values); 
    }

    // Retrieve region_boundary
    std::vector<int> region_boundary;
    for (const auto& boundary_index : pt.get_child("region_boundary")) {
        region_boundary.push_back(boundary_index.second.get_value<int>());
    }

    // Create inputData struct and populate it
    InputData input_data;
    input_data.points = points;
    input_data.region_boundary = region_boundary;
    input_data.additional_constraints = additional_constraints;

    return input_data; 
}
