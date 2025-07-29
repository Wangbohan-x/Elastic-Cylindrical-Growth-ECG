#include"grid.h"
#include"m_cylinder.h"
#include"polygon.h"

bool OptimizeCylinder(const std::vector<Point>& points,
	Vec3f& axis_dir_in,
	std::vector<Vec3f>& axis_point_in,
	double& radius,
	int tag
	/*double& optimized_s0,
	double& optimized_s1*/);

//写为ply
void write_ply(Point_set& points, std::string path) {
	std::ofstream out(path);
	CGAL::IO::set_ascii_mode(out);
	out << points;
	//每个类别记为一个label  没有的设为0
	return;
}

bool is_direction_near_x_axis(const Vec3f& dir, double max_angle_deg) {
	Vec3f x_axis_pos(1.0, 0.0, 0.0);
	Vec3f x_axis_neg(-1.0, 0.0, 0.0);

	double dot_pos = dir.dot(x_axis_pos);
	double dot_neg = dir.dot(x_axis_neg);
	double angle_rad = std::acos(max(std::abs(dot_pos), std::abs(dot_neg)));
	//std::cout << angle_rad*180/ 3.141592654 <<" ";
	double max_angle_rad = max_angle_deg * 3.141592654 / 180.0;
	return angle_rad < max_angle_rad;
}

bool is_direction_near_z_axis(const Vec3f& dir, double max_angle_deg) {
	Vec3f z_axis_pos(0.0, 0.0, 1.0);
	Vec3f z_axis_neg(0.0, 0.0, -1.0);

	double dot_pos = dir.dot(z_axis_pos);
	double dot_neg = dir.dot(z_axis_neg);
	double angle_rad = std::acos(max(std::abs(dot_pos), std::abs(dot_neg)));
	//std::cout << angle_rad*180/ 3.141592654 <<" ";
	double max_angle_rad = max_angle_deg * 3.141592654 / 180.0;
	return angle_rad < max_angle_rad;
}

double calculate_distance(const Vec3f& p1, const Vec3f& p2) {
	double dx = p2.X() - p1.X();
	double dy = p2.Y() - p1.Y();
	double dz = p2.Z() - p1.Z();
	double ret = sqrt(dx * dx + dy * dy + dz * dz);
	//std::cout << ret<<"   ";
	return ret;
}

// 添加颜色属性
void assign_color_map(Point_set& points) {
	Point_set::Property_map<unsigned char>
		red = points.add_property_map<unsigned char>("red", 255).first,
		green = points.add_property_map<unsigned char>("green", 255).first,
		blue = points.add_property_map<unsigned char>("blue", 255).first;
	return;
}


//添加某种颜色
void assign_uniform_color(Point_set& points) {
	auto red = points.property_map<unsigned char>("red").first;
	auto green = points.property_map<unsigned char>("green").first;
	auto blue = points.property_map<unsigned char>("blue").first;

	const unsigned char r = static_cast<unsigned char>(rand() % 255);
	const unsigned char g = static_cast<unsigned char>(rand() % 255);
	const unsigned char b = static_cast<unsigned char>(rand() % 255);
	for (auto id : points) {
		put(red, id, r);
		put(green, id, g);
		put(blue, id, b);
	}
	return;
}

// 如果不存在 label 属性，则添加
void add_label_property(Point_set& points) {
	
	if (!points.property_map<int>("label").second) {
		points.add_property_map<int>("label", -1);
	}
}

void set_label_for_indices(Point_set& points,
	const Region_growing::Region indices,
	int label_value) 
{
	// 确保 label 属性已存在
	auto label = points.property_map<int>("label").first;
	if (!points.property_map<int>("label").second) {
		std::cerr << "Error: 'label' property not found. Please call add_label_property() first." << std::endl;
		return;
	}
	
	for (auto& id:indices) {
		label[id] = label_value;
	}
	return;
}

void set_label(Point_set& points,
	int label_value)
{
	auto label = points.property_map<int>("label").first;
	if (!points.property_map<int>("label").second) {
		std::cerr << "Error: 'label' property not found. Please call add_label_property() first." << std::endl;
		return;
	}

	for (auto& id : points) {
		label[id] = label_value;
	}
	return;
}


// 合并点集
void merge_point(Point_set& target, const Point_set& source) {
	for (auto idx : source) {
		target.insert(source.point(idx));
	}
}

std::vector<Point> ExtractPointsFromPointSet(const Point_set& point_set) {
	std::vector<Point> points;
	points.reserve(point_set.size());

	for (auto it = point_set.begin(); it != point_set.end(); ++it) {
		points.push_back(point_set.point(*it));
	}

	return points;
}

std::vector<Point> ConvertRegionToPoints(
	const Region_growing::Region& region,
	const Point_set& point_set)
{
	std::vector<Point> result;
	result.reserve(region.size());

	for (std::size_t idx : region) {
		result.push_back(point_set.point(idx));
	}

	return result;
}

// 检测圆柱体1.0
std::pair<std::vector<CylinderElement>, std::vector<Point_set>>
detect_cylinders(const Point_set& input,
	int k = 20,
	double max_distance = 0.008,//0.02
	double max_angle = 90.0,
	std::size_t min_region_size = 20) {//60
	std::vector<Point_set>point_sets;
	std::vector<CylinderElement>cylinders;

	assert(input.has_normal_map());
	Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(input, CGAL::parameters::k_neighbors(k));
	Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_cylinder_fit_region(
		input,
		CGAL::parameters::
		maximum_distance(max_distance).
		maximum_angle(max_angle).
		minimum_region_size(min_region_size));
	Region_growing region_growing(input, neighbor_query, region_type);
	std::vector<typename Region_growing::Primitive_and_region> regions;
	region_growing.detect(boost::make_function_output_iterator(
		[&](const Primitive_and_region& region) {
			Point_set temPnts;
			for (auto id : region.second) {
				temPnts.insert(input, id);
			}
			point_sets.emplace_back(temPnts);
			CylinderElement cyl;

			const auto& cylinder_param = region.first;
			const auto& dx = cylinder_param.axis.direction().dx();
			const auto& dy = cylinder_param.axis.direction().dy();
			const auto& dz = cylinder_param.axis.direction().dz();
			const auto& cx = cylinder_param.axis.point(0).x();
			const auto& cy = cylinder_param.axis.point(0).y();
			const auto& cz = cylinder_param.axis.point(0).z();
			const auto& ra = cylinder_param.radius;
			double height = GetHeight(temPnts, dx, dy, dz, cx, cy, cz);

			cyl.radius = ra;
			cyl.direction = Vec3f(dx, dy, dz);
			Vec3f position(cx, cy, cz);
			Vec3f position1 = position + 0.5 * height * cyl.direction;
			Vec3f position2 = position - 0.5 * height * cyl.direction;
			cyl.positions.push_back(position);//中心点  下面是两个端点
			cyl.positions.push_back(position1);
			cyl.positions.push_back(position2);
			cylinders.push_back(cyl);
		}
	)
	);
	return { cylinders, point_sets };
}

// 检测圆柱体2.0
std::pair<std::vector<CylinderElement>, std::vector<Region_growing::Region>>
detect_cylinders2(const Point_set& input,
	int k = 20,
	double max_distance = 0.007,//0.02
	double max_angle = 90.0,
	std::size_t min_region_size = 10) {//60
	//std::vector<Point_set>point_sets;
	std::vector<Region_growing::Region>regionss;
	std::vector<CylinderElement>cylinders;

	assert(input.has_normal_map());
	Neighbor_query neighbor_query = CGAL::Shape_detection::Point_set::make_k_neighbor_query(input, CGAL::parameters::k_neighbors(k));
	Region_type region_type = CGAL::Shape_detection::Point_set::make_least_squares_cylinder_fit_region(
		input,
		CGAL::parameters::
		maximum_distance(max_distance).
		maximum_angle(max_angle).
		minimum_region_size(min_region_size));
	Region_growing region_growing(input, neighbor_query, region_type);
	std::vector<typename Region_growing::Primitive_and_region> regions;
	region_growing.detect(boost::make_function_output_iterator(
		[&](const Primitive_and_region& region) {
			Point_set temPnts;
			for (auto id : region.second) {
				temPnts.insert(input, id);
			}
			//point_sets.emplace_back(temPnts);
			regionss.emplace_back(region.second);
			CylinderElement cyl;

			const auto& cylinder_param = region.first;
			const auto& dx = cylinder_param.axis.direction().dx();
			const auto& dy = cylinder_param.axis.direction().dy();
			const auto& dz = cylinder_param.axis.direction().dz();
			const auto& cx = cylinder_param.axis.point(0).x();
			const auto& cy = cylinder_param.axis.point(0).y();
			const auto& cz = cylinder_param.axis.point(0).z();
			const auto& ra = cylinder_param.radius;
			double height = GetHeight(temPnts, dx, dy, dz, cx, cy, cz);

			cyl.radius = ra;
			cyl.direction = Vec3f(dx, dy, dz);
			Vec3f position(cx, cy, cz);
			Vec3f position1 = position + 0.5 * height * cyl.direction;
			Vec3f position2 = position - 0.5 * height * cyl.direction;
			cyl.positions.push_back(position);//中心点  下面是两个端点
			cyl.positions.push_back(position1);
			cyl.positions.push_back(position2);
			cylinders.push_back(cyl);
		}
	)
	);
	//return { cylinders, point_sets };
	return{ cylinders,regionss };
}

//球体区域的点
Point_set collect_points_in_sphere(
	const Point_set& all_points,
	const CGAL::Point_3<Kernel>& center,
	double radius)
{
	Point_set result;
	double sqr_radius = radius * radius;
	for (auto id : all_points) {
		const auto& p = all_points.point(id);
		if (CGAL::squared_distance(p, center) <= sqr_radius) {
			result.insert(p);
		}
	}
	return result;
}

Vec3f toVec3f(const CGAL::Point_3<Kernel>& pt) {
	return Vec3f(pt.x(), pt.y(), pt.z());
}

//圆柱体范围的点
Point_set sample_cylinder_region(
	const Point_set& all_points,
	const Vec3f& endpoint,
	const Vec3f& axis_dir,
	Region_growing::Region& ids,
	double radius,
	double height)
{
	Point_set sampled;

	//Vec3f axis_dir(1, 0, 0); // x 轴正方向，已归一化
	//Vec3f axis_dir(0, 0, -1); // z 轴负方向
	Vec3f base = endpoint;
	Vec3f top = endpoint + height * axis_dir;

	for (auto idx : all_points) {
		Vec3f pt = toVec3f(all_points.point(idx));

		// 投影到轴方向上的距离
		double proj_len = (pt - base).dot(axis_dir);

		// 判断是否在圆柱高度范围内
		if (proj_len < 0 || proj_len > height) continue;

		// 计算点到轴的垂直距离
		Vec3f proj_point = base + proj_len * axis_dir;
		double radial_dist = (pt - proj_point).length();

		if (radial_dist <= radius) {
			sampled.insert(all_points.point(idx), all_points.normal(idx));
			ids.emplace_back(idx);
		}
	}

	return sampled;
}

bool grow_cylinder_one_end(
	CylinderElement& initial_cylinder,
	const Point_set& all_points,
	/*Point_set&ret,*/
	Region_growing::Region&ret,
	int num,
	int outcount,
	double grow_step,
	double max_distance,
	double max_angle,
	double max_radius_diff
	){
	Vec3f direction = initial_cylinder.direction;
	Vec3f head = initial_cylinder.positions[1];
	Vec3f tail = initial_cylinder.positions[2];

	//head = head.X() > tail.X() ? head : tail;
	//direction = direction.X() > 0 ? direction : -direction;
	//Vec3f axis(1, 0, 0);

	head = head.Z() < tail.Z() ? head : tail;
	direction = direction.Z() < 0 ? direction : -direction;
	Vec3f axis(0, 0, -1);
	double dis = 0.008;
	if (num>0)
	{
		direction = axis;
		dis = dis + num*0.001;
	}
	//std::cout <<  direction.X() << ", " << direction.Y()<< ", " << direction.Z() << "\n";
	
	//head = head.Z() < tail.Z() ? head : tail;
	//dir.Z() = dir.Z() < 0 ? dir.Z() : -dir.Z();
	//Vec3f axis(0, 0, -1);
	//Vec3f new_center = head + dir * grow_step;
	//CGAL::Point_3<Kernel> center(new_center[0], new_center[1], new_center[2]);

	//Point_set region = collect_points_in_sphere(all_points, center, grow_step);// 收集球体内点
	Region_growing::Region ids;
	Point_set region = sample_cylinder_region(all_points, head+num* direction */*grow_step*/0.03, direction, ids, 0.03+num*0.005, grow_step);// 收集圆柱体内点

	//write_ply(region, "C:/Users/Shawn/Desktop/CylinderDetection/RM/region.ply");
	
	// 如果没有检测到圆柱体，就缩小圆柱体范围采样 然后对其中的点进行圆柱体的拟合
	//if (num > 0)
	//{
	//	std::cout << "   " << num<< "    ";
	//	write_ply(region, "C:/Users/Shawn/Desktop/CylinderDetection/RC/region.ply");
	//}
	write_ply(region, "C:/Users/Shawn/Desktop/CylinderDetection/RM/region.ply");
	if (region.size() < 5) 
	{
		return false;
	}
	
	region.add_normal_map();
	CGAL::jet_estimate_normals<CGAL::Sequential_tag>(region, 12, region.parameters().degree_fitting(2));

	//auto candidates = detect_cylinders(region);
	auto candidates = detect_cylinders2(region, 20, dis, 90, 10);
	//int i = 0;
	//for (auto points:candidates.second)
	//{
	//	assign_uniform_color(points);
	//	std::string pathtemp = "C:/Users/Shawn/Desktop/CylinderDetection/RC/"+std::to_string(i)+".ply";
	//	std::ofstream out(pathtemp);
	//	CGAL::IO::set_ascii_mode(out);
	//	out << points;
	//	i++;
	//}
	if (candidates.first.size()>1)
	{
		//按照点多少排序
		std::vector<std::pair<CylinderElement, Region_growing::Region>> combined;
		for (std::size_t i = 0; i < candidates.first.size(); ++i) {
			combined.emplace_back(candidates.first[i], candidates.second[i]);
		}

		// 按照 region 点数量排序（降序）
		std::sort(combined.begin(), combined.end(),
			[](const auto& a, const auto& b) {
				return a.second.size() > b.second.size();  // 根据 Region 点数降序排列
			});

		// 解包回原来的两个 vector
		for (std::size_t i = 0; i < combined.size(); ++i) {
			candidates.first[i] = combined[i].first;
			candidates.second[i] = combined[i].second;
		}
	}
	for (int i = 0; i < candidates.first.size(); ++i) {
		//std::cout << candidates.first[i].radius	 << std::endl;
		auto cylpoint = ConvertRegionToPoints(candidates.second[i], all_points);
		OptimizeCylinder(cylpoint, candidates.first[i].direction, candidates.first[i].positions, candidates.first[i].radius,1);
		auto& prim = candidates.first[i];
		Vec3f dir(prim.direction.X(), prim.direction.Y(), prim.direction.Z());
		Vec3f position1 = prim.positions[1];
		Vec3f position2 = prim.positions[2];
		double r = prim.radius;
		double h = calculate_distance(position1, position2);
		//std::cout <<r<<" " <<h<< std::endl;
		position1 = position1.X() > position2.X() ? position2 : position1;

		//position1 = position1.Z() > position2.Z() ? position1 : position2;
		//if (!is_direction_near_x_axis(dir, 15)) {
		//	continue;
		//}

		if (!is_direction_near_z_axis(dir, 20)) {
			continue;
		}
		//std::cout << calculate_distance(position1, head) << std::endl;
		if (calculate_distance (position1,head)<0.2+num*grow_step)//0.2
		{
			//std::cout << 3 << " ";
			//merge_point_sets_simple(ret, candidates.second[i]);
			//ret.insert(ret.end(),candidates.second[i].begin(), candidates.second[i].end());
			Point_set temp;

			for (auto j : candidates.second[i])
			{
				ret.emplace_back(ids[j]);
				temp.insert(region.point(j));
			}
			//std::cout << candidates.first[i].radius << std::endl;
			//assign_color_map(temp);
			//assign_uniform_color(temp);
			//write_ply(temp,"C:/Users/Shawn/Desktop/CylinderDetection/RC/ret/test" + std::to_string(outcount) + ".ply");


			//std::cout << ret.size()<<std::endl;
			initial_cylinder = prim;
			return true;
		}
	}
	//current_cylinder.positions.push_back(endpoint);
	std::cout << "没找到合适的圆柱体" << std::endl;
	return false;
}

//弹性生长
/*Point_set*/Region_growing::Region  grow_cylinder_elastically(
	CylinderElement& initial_cylinder,
	const Point_set& all_points,
	int max_iters = 100 ,
	double grow_step = 0.15,
	double max_distance = 0.02,
	double max_angle = 10.0,
	double max_radius_diff = 0.005) {
	//double radius = initial_cylinder.radius;
	
	//Point_set ret;
	Region_growing::Region ret;
	int num = 0;
	for (int i = 0; i < max_iters; ++i) {

		bool head_grown = grow_cylinder_one_end(initial_cylinder, all_points, ret ,num,i,grow_step, max_distance, max_angle, max_radius_diff);
		//std::cout << num;
		//bool tail_grown = grow_cylinder_one_end(tail, -dir, all_points, grow_radius, max_distance, max_angle, max_radius_diff, grown_cylinder);
		if (head_grown)
		{
			num = 0;
		}
		if (!head_grown/* && !tail_grown*/) {

			num++;
		}
		if (num > 10)break;

	}

	return ret;

}

std::pair<CylinderElement, Point_set>
select_largest_cylinder(const std::pair<std::vector<CylinderElement>, std::vector<Point_set>>& detected) {
	const auto& cylinders = detected.first;
	const auto& point_sets = detected.second;

	assert(cylinders.size() == point_sets.size());
	if (cylinders.empty()) {
		throw std::runtime_error("No cylinders detected.");
	}

	std::size_t max_index = 0;
	std::size_t max_size = point_sets[0].size();

	for (std::size_t i = 1; i < point_sets.size(); ++i) {
		if (point_sets[i].size() > max_size) {
			max_index = i;
			max_size = point_sets[i].size();
		}
	}

	return { cylinders[max_index], point_sets[max_index] };
}



void export_cylinder_axes_to_txt(std::vector<CylinderElement>& cylinders, const std::string& filename) {
	std::ofstream ofs(filename);
	if (!ofs.is_open()) {
		throw std::runtime_error("Cannot open output file.");
	}

	ofs << std::fixed << std::setprecision(6);  // 控制浮点输出精度

	for (auto& cyl : cylinders) {
		if (cyl.positions.size() < 2) continue;  

		const Vec3f& origin = cyl.positions[0];       
		Vec3f& direction = cyl.direction;       // 单位方向向量
		if (direction.X() < 0) {
			direction = -direction;
		}
		ofs << origin.X() << " " << origin.Y()<< " " << origin.Z()<< " "
			<< direction.X() << " " << direction.Y() << " " << direction.Z() << "\n";
	}

	ofs.close();
}

void save_radius_to_txt(const std::string& file_path, double radius) {
	std::ofstream out(file_path);
	if (!out.is_open()) {
		std::cerr << "Failed to open file: " << file_path << std::endl;
		return;
	}
	out << radius << std::endl;
	out.close();
}

void compute_weighted_average_radius(const std::string& file_path,const std::vector<CylinderElement>& cylinders) {
	double weighted_sum = 0.0;
	double total_length = 0.0;

	for (const auto& cyl : cylinders) {
		assert(cyl.positions.size() >= 2);
		double length = calculate_distance(cyl.positions[1], cyl.positions[2]);
		weighted_sum += cyl.radius * length;
		total_length += length;
	}
	if (total_length == 0.0) return;
	double ret = weighted_sum / total_length;
	//std::cout << std::endl;
	//std::cout << ret << std::endl;
	save_radius_to_txt(file_path, ret);
	return;
}

struct RadialResidual {
	RadialResidual(const Point& p) : p_(p) {}

	template <typename T>
	bool operator()(const T* const axis_dir,      // 单位方向向量 (3)
		const T* const axis_point,    // 轴线上一点 (3)
		const T* const radius,        // 半径 (1)
		T* residual) const {
		// 计算向量 P - P0
		T wx = T(p_.x()) - axis_point[0];
		T wy = T(p_.y()) - axis_point[1];
		T wz = T(p_.z()) - axis_point[2];
		// 计算投影长度 t
		T t = wx * axis_dir[0] + wy * axis_dir[1] + wz * axis_dir[2];
		// 最近点 C
		T cx = axis_point[0] + t * axis_dir[0];
		T cy = axis_point[1] + t * axis_dir[1];
		T cz = axis_point[2] + t * axis_dir[2];
		// 径向距离
		T dx = T(p_.x()) - cx;
		T dy = T(p_.y()) - cy;
		T dz = T(p_.z()) - cz;
		T dist = sqrt(dx * dx + dy * dy + dz * dz);
		residual[0] = dist - radius[0];
		return true;
	}

private:
	const Point p_;
};

double ComputeAngleBetweenDirections( Vec3f& v1,  Vec3f& v2) {
	double dot = v1.dot(v2);  // 假设 * 是点乘运算符
	double norm1 = v1.normalize();
	double norm2 = v2.normalize();

	if (norm1 == 0.0 || norm2 == 0.0) {
		// 如果任一向量长度为 0，夹角无意义，返回一个非法值
		return -1.0;
	}

	double cos_theta = dot / (norm1 * norm2);

	// 限制 cos_theta 的范围在 [-1, 1] 之间，避免 acos 出现非法输入
	cos_theta = std::max(-1.0, std::min(1.0, cos_theta));
	double angle_deg = std::acos(cos_theta) * 180.0 / 3.14159265;
	if (angle_deg > 90.0) angle_deg = 180.0 - angle_deg;
	// 转为角度
	return angle_deg;
}


bool OptimizeCylinder(const std::vector<Point>& points,
	Vec3f& axis_dir_in,
	std::vector<Vec3f>& axis_point_in,
	double& radius,
	int tag
	/*double& optimized_s0,
	double& optimized_s1*/) {
	double axis_dir[3] = { axis_dir_in.X(), axis_dir_in.Y(), axis_dir_in.Z() };
	double axis_point[3] = { axis_point_in[0].X(), axis_point_in[0].Y(), axis_point_in[0].Z() };
	double r = radius;
	// Ceres 问题定义
	ceres::Problem problem;
	// 侧面残差
	for (const auto& p : points) {
		problem.AddResidualBlock(
			new ceres::AutoDiffCostFunction<RadialResidual, 1, 3, 3, 1>(
				new RadialResidual(p)),
			nullptr,
			axis_dir, axis_point, &radius);
	}
	
	ceres::LocalParameterization* unit_vec =
		new ceres::HomogeneousVectorParameterization(3);
	problem.SetParameterization(axis_dir, unit_vec);

	// Solver 参数
	ceres::Solver::Options options;
	options.linear_solver_type = ceres::DENSE_QR;
	options.minimizer_progress_to_stdout = false;
	ceres::Solver::Summary summary;
	ceres::Solve(options, &problem, &summary);

	if (!summary.IsSolutionUsable()) {
		std::cerr << "Optimization failed: " << summary.BriefReport() << std::endl;
		return false;
	}

	auto v = Vec3f(axis_dir[0], axis_dir[1], axis_dir[2]);
	double theta = ComputeAngleBetweenDirections(v, axis_dir_in);
	//std::cout << theta << std::endl;
	double min_t = std::numeric_limits<double>::infinity();
	double max_t = -std::numeric_limits<double>::infinity();
	for (const auto& p : points) {
		double wx = p.x() - axis_point[0];
		double wy = p.y() - axis_point[1];
		double wz = p.z() - axis_point[2];
		double t = wx * axis_dir[0] + wy * axis_dir[1] + wz * axis_dir[2];
		if (t < min_t) min_t = t;
		if (t > max_t) max_t = t;
	}
	//optimized_s0 = min_t;
	//optimized_s1 = max_t;
	if (fabs(radius - r) > 0.2 * r)
	{
		radius = r;
	}
	//if (theta>20) {
	//	return true;
	//}
	axis_dir_in = Vec3f(axis_dir[0], axis_dir[1], axis_dir[2]);
	if (tag == 1)return true;
	axis_point_in[0] = Vec3f(axis_point[0], axis_point[1], axis_point[2]);
	axis_point_in[1] = axis_point_in[0] + min_t * axis_dir_in;
	axis_point_in[2] = axis_point_in[0] + max_t * axis_dir_in;


	return true;
}

void Optimize(CylinderElement &cyl , const Point_set points) {

	OptimizeCylinder(ExtractPointsFromPointSet(points), cyl.direction, cyl.positions, cyl.radius,0);
	return;
}

int main(int argc, char** argv){


	std::string path = "C:/Users/Shawn/Desktop/CylinderDetection/RM/";
	Point_set base;
	Point_set other;
	CGAL::IO::read_LAS(path+"3mmZbase.las" , base);
	CGAL::IO::read_LAS(path + "3mmother.las", other);
	std::cout << base.size() << std::endl;
	std::cout << other.size() << std::endl;


	base.add_normal_map();
	other.add_normal_map();

	assert(base.has_normal_map());
	assert(other.has_normal_map());

	CGAL::jet_estimate_normals<CGAL::Sequential_tag>(base, 12, base.parameters().degree_fitting(2));  
	CGAL::jet_estimate_normals<CGAL::Sequential_tag>(other, 12, other.parameters().degree_fitting(2));  



		std::vector<Point>points;
		std::vector<std::vector<size_t>>polys;

		std::vector<Point>points1;
		std::vector<std::vector<size_t>>polys1;

		add_label_property(base);
		add_label_property(other);
		auto base_ret = detect_cylinders(base,20,0.02,90,500);          
		std::cout << base_ret.second.size();	


		for (int i = 0; i<base_ret.first.size();i++)
		{		
			std::cout <<std::endl << i << std::endl;

			//set_label(base_ret.second[i], i);
			//write_ply(base, "C:/Users/Shawn/Desktop/CylinderDetection/RC/test.ply");
			auto ininf = base_ret.first[i];
			auto basepoint = base_ret.second[i];
			Optimize(ininf, basepoint);

			Region_growing::Region other_ret = grow_cylinder_elastically(ininf, other);
			std::cout << std::endl << other_ret.size() << std::endl;
			set_label_for_indices(other, other_ret,i);


			
			
			Point_set other_out;
			for (auto id: other_ret)
			{
				other_out.insert(other.point(id));
			}
			assign_color_map(other_out);
			assign_uniform_color(other_out);
			write_ply(other_out, path+"retz/other_out"+ std::to_string(i) +".ply");
			
			
			int id_begin = other.size();
			merge_point(other, base_ret.second[i]);
			auto label = other.property_map<int>("label").first;
			for (size_t id = id_begin,j = 0; j < base_ret.second[i].size(); id++,j++)
			{
				label[id] = i;
			}
			//assign_uniform_color(other_ret);
			//write_ply(other_ret, "C:/Users/Shawn/Desktop/CylinderDetection/RC/other_ret"+std::to_string(i)+".ply");

		}

		write_ply(other, path+"retz/ret.ply");

	return 0;
}