#ifndef UTILS_H
#define UTILS_H

#include <CGAL/Delaunay_triangulation_2.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>

typedef CGAL::Exact_predicates_exact_constructions_kernel K;
typedef CGAL::Constrained_Delaunay_triangulation_2<K> DT;
typedef DT::Face_handle FaceHandle;
typedef DT::Point Point;

template <typename P>
double angle_between(const P& p1, const P& p2, const P& p3);

template <typename FaceHandle>
int obtuse_vertex_index(const FaceHandle& face);

int count_obtuse_triangles(const DT& dt);

template <typename DT>
std::vector<std::pair<size_t, size_t>> print_edges(const DT& dt, std::vector<Point> points);

#endif // UTILS_H
