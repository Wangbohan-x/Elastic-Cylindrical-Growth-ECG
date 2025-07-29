#ifndef SELF_UTILS_H
#define SELF_UTILS_H

#define CGAL_LINKED_WITH_LASLIB
#define GLOG_USE_GLOG_EXPORT
#define GLOG_NO_ABBREVIATED_SEVERITIES
#include <fstream>
#include <limits>
#include<algorithm>
#include <utility>
#include <vector>

#include <CGAL/Point_set_3.h>
#include<CGAL/Shape_detection/Efficient_RANSAC.h>
#include<CGAL/Shape_detection/Region_growing.h> //这个要使用eigen库
#include <CGAL/Point_set_3/IO/LAS.h>
#include <CGAL/Exact_predicates_inexact_constructions_kernel.h>
#include <CGAL/jet_estimate_normals.h>
#include <CGAL/IO/PLY.h>

#include<CGAL/Simple_cartesian.h>

#include<glog/logging.h>
#include<ceres/ceres.h>
#include <ceres/rotation.h>

//typedef CGAL::Exact_predicates_inexact_constructions_kernel Kernel;
//using Kernel = CGAL::Exact_predicates_inexact_constructions_kernel;
using Kernel = CGAL::Simple_cartesian<double>;
//typedef Kernel::Point_3 Point;
using  Point = Kernel::Point_3;
using  Vector = Kernel::Vector_3;

//typedef CGAL::Point_set_3<Point> Point_set;
using Point_set = CGAL::Point_set_3<Point,Vector>;
using FT = Kernel::FT;
using Neighbor_query = CGAL::Shape_detection::Point_set::K_neighbor_query_for_point_set<Point_set>;
using Region_type = CGAL::Shape_detection::Point_set::Least_squares_cylinder_fit_region_for_point_set<Point_set>;
using Region_growing = CGAL::Shape_detection::Region_growing<Neighbor_query, Region_type>;
typedef typename Region_growing::Primitive_and_region Primitive_and_region;


#endif 