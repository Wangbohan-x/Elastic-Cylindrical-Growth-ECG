#pragma once
#include"utils.h"
#include"m_cylinder.h"
void CylinderPolygon(CylinderElement& cylinder,std::vector<std::vector<size_t>>&polys,std::vector<Point>&points,int circle_point_num = 10) {
	double pi = 4.0 * atan(1.0);
	double radius = cylinder.radius;
	
	Vec3f v = Vec3f(1.0f, 0.0f, 0.0f);
	Vec3f u = v.cross(cylinder.direction);
	v = u.cross(cylinder.direction);
	u.normalize();
	v.normalize();
	int num = points.size();
	int  begin;
	std::vector<size_t>axis;
	points.push_back(Point(cylinder.positions[0].X(), cylinder.positions[0].Y(), cylinder.positions[0].Z()));
	axis.push_back(num);
	num++;
	for (size_t i = 1; i < cylinder.positions.size(); i++)
	{
		begin = num;
		std::vector<size_t>tmp;
		Vec3f Position = cylinder.positions[i];

		//std::cout << Position.X() << " "<< Position.Y()<< "  " << Position.Z()<<"   ";

		double phi_inc = 2 * pi / circle_point_num;
		double phi = 0.0;
		points.push_back(Point(Position.X(), Position.Y(), Position.Z()));
		axis.push_back(num);
		num++;
		
		for (int j = 0; j < circle_point_num; ++j, ++num) {
			tmp.push_back(num);
			Vec3f Tempoint = Position + radius * (cos(phi) * u + sin(phi) * v);
			//std::cout << Tempoint << "\t";
			points.push_back(Point(Tempoint.X(), Tempoint.Y(), Tempoint.Z()));
			phi += phi_inc;
		}
		polys.push_back(tmp);
		


		//for (size_t i = 0; i < circle_point_num; i++)
		//{
		//	std::vector<size_t>tri;
		//	tri.push_back(begin);
		//	if (i == circle_point_num-1)
		//	{
		//		tri.push_back(tmp[i]);
		//		tri.push_back(tmp[0]);
		//	}
		//	else {
		//		tri.push_back(tmp[i]);
		//		tri.push_back(tmp[i + 1]);
		//	}
		//	polys.push_back(tri);
		//}
	}
	polys.push_back(axis);

	//std::vector<size_t>tmp;
	//for (size_t i = 0; i < cylinder.positions.size(); i++)
	//{
	//	points.push_back(Point(cylinder.positions[i].X(), cylinder.positions[i].Y(), cylinder.positions[i].Z()));
	//	tmp.push_back(num);
	//	num++;
	//}
	//polys.push_back(tmp);
	
	
}