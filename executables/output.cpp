#define BOOST_BIND_GLOBAL_PLACEHOLDERS
#include <boost/bind.hpp>
#include <boost/property_tree/ptree.hpp>
#include <boost/property_tree/json_parser.hpp>
#include <iostream>
#include <vector>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include "output.h"
#include "inputs.h"
#include <string>
#include <sstream> 
#include <boost/algorithm/string/replace.hpp>
#include <fstream>

// Define CGAL types
typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef K::Point_2 Point;

void write_json_no_escaping(const boost::property_tree::ptree& pt, const std::string& filename) {
    std::ostringstream oss;
    write_json(oss, pt, true); 

    std::string json_str = oss.str();
    
    // Replace the escaped `\/` with `/`
    boost::replace_all(json_str, "\\/", "/");

    // Write the modified JSON string to the file
    std::ofstream outfile(filename);
    if (outfile) {
        outfile << json_str;
        outfile.close();
        std::cout << "Output written to " << filename << std::endl;
    } else {
        std::cerr << "Error opening file: " << filename << std::endl;
    }
}

// Function to convert K::FT to a string representation
std::string rational_to_string(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    if (exact_coord.get_den() == 1) {
        return exact_coord.get_num().get_str();
    } else {
        std::ostringstream oss;
        oss << exact_coord.get_num().get_str() << "/" << exact_coord.get_den().get_str();
        return oss.str();
    }
}

void print_rational(const K::FT& coord) {
    const auto exact_coord = CGAL::exact(coord);
    std::cout << exact_coord.get_num() << "/" << exact_coord.get_den();
}

void output(const std::vector<std::pair<Point, Point>>& edges, std::vector<Point> steiner_points_given) {
    // Creation of property tree
    boost::property_tree::ptree pt;

    // Read from JSON file
    try {
        read_json("../input.json", pt);
    } catch (const boost::property_tree::json_parser_error &e) {
        std::cerr << "Error reading JSON: " << e.what() << std::endl;
    }

    // Get data from property tree
    std::string instance_uid = pt.get<std::string>("instance_uid");

    std::vector<Point> steiner_points = steiner_points_given;

    // Create the output property tree
    boost::property_tree::ptree output_pt;

    // Populate the JSON structure
    output_pt.put("content_type", "CG_SHOP_2025_Solution");
    output_pt.put("instance_uid", instance_uid);

    // Steiner points x
    boost::property_tree::ptree steiner_points_x_node;
    for (const auto& x : steiner_points) {
        boost::property_tree::ptree child;
        child.put("", rational_to_string(x.x()));
        steiner_points_x_node.push_back(std::make_pair("", child));
    }
    output_pt.add_child("steiner_points_x", steiner_points_x_node);

    // Steiner points y
    boost::property_tree::ptree steiner_points_y_node;
    for (const auto& y : steiner_points) {
        boost::property_tree::ptree child;
        child.put("", rational_to_string(y.y()));
        steiner_points_y_node.push_back(std::make_pair("", child));
    }
    output_pt.add_child("steiner_points_y", steiner_points_y_node);

    // Edges
    boost::property_tree::ptree edges_node;
    for (const auto& edge : edges) {
        boost::property_tree::ptree edge_node;
        
        // For each edge's first point
        boost::property_tree::ptree edge_first;
        edge_first.put("x", rational_to_string(edge.first.x()));
        edge_first.put("y", rational_to_string(edge.first.y()));

        // For each edge's second point
        boost::property_tree::ptree edge_second;
        edge_second.put("x", rational_to_string(edge.second.x()));
        edge_second.put("y", rational_to_string(edge.second.y()));

        // Add first and second points to the edge node
        edge_node.add_child("first", edge_first);
        edge_node.add_child("second", edge_second);

        // Add the edge node to edges_node
        edges_node.push_back(std::make_pair("", edge_node));
    }
    // Add edges_node to the main output property tree
    output_pt.add_child("edges", edges_node);

    // Write the output JSON to a file
    try {
        write_json_no_escaping(output_pt, "../output.json");
    } catch (const boost::property_tree::json_parser_error &e) {
        std::cerr << "Error writing JSON: " << e.what() << std::endl;
    }
}
