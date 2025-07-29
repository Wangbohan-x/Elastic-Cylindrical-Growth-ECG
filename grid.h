#pragma once
#include"utils.h"

#include<CGAL/Bbox_3.h>

class Grid
{
public:
	Grid(Point_set& point_set,double step):step(step){
		RecalculateAABB(point_set);
		this->splitX = ceil((this->aabb.x_span()) / step);
		this->splitY = ceil((this->aabb.y_span()) / step);
		this->splitZ = ceil((this->aabb.z_span()) / step);
		this->aabb.rep[3] += splitX * step - this->aabb.x_span();
		this->aabb.rep[4] += splitY * step - this->aabb.y_span();
		this->aabb.rep[5] += splitZ * step - this->aabb.z_span();
		std::cout << "Grid splitX: " << splitX << std::endl;
		std::cout << "Grid splitY: " << splitY << std::endl;
		std::cout << "Grid splitZ: " << splitZ << std::endl;
		std::cout << "AABB:" << std::endl << this->aabb << std::endl;
		BuildGrid(point_set);
		std::cout << "Grid num: " << idx_data.size() << std::endl;
	}
	~Grid() {
		
	}

private:
	void RecalculateAABB(Point_set& point_set) {
		double x, y, z;
		x = point_set.point(0).x();
		y = point_set.point(0).y();
		z = point_set.point(0).z();
		this->aabb.rep[0] = x;
		this->aabb.rep[1] = y;
		this->aabb.rep[2] = z;
		this->aabb.rep[3] = x;
		this->aabb.rep[4] = y;
		this->aabb.rep[5] = z;
		for (size_t i = 1; i < point_set.size(); i++)
		{
			x = point_set.point(i).x();
			y = point_set.point(i).y();
			z = point_set.point(i).z();
			this->aabb.rep[0] = min(this->aabb.rep[0], x);
			this->aabb.rep[1] = min(this->aabb.rep[1], y);
			this->aabb.rep[2] = min(this->aabb.rep[2], z);
			this->aabb.rep[3] = max(this->aabb.rep[3], x);
			this->aabb.rep[4] = max(this->aabb.rep[4], y);
			this->aabb.rep[5] = max(this->aabb.rep[5], z);
		}
	}

	void BuildGrid(Point_set& point_set) {
		long long x, y, z;
		long long key;
		std::vector<size_t> v;
		v.clear();
		for (size_t i = 0; i < point_set.size(); i++) {
			x = floor((point_set.point(i).x() - aabb.xmin()) / this->step);
			y = floor((point_set.point(i).y() - aabb.ymin()) / this->step);
			z = floor((point_set.point(i).z() - aabb.zmin()) / this->step);
			key = x << 40 | y << 20 | z;
			if (idx_data.find(key) != idx_data.end())idx_data[key].emplace_back(i);
			else idx_data.emplace(key, v);
		}
	}



public:
	void GetAllGridPoints(std::vector<std::vector<size_t>>& indices) {
		for (auto it = idx_data.begin(); it != idx_data.end(); it++)
		{
			indices.emplace_back(it->second);
		}
	}

//private:

	CGAL::Bbox_3 aabb;
	double step;
	int splitX, splitY, splitZ;
	std::map<long long, std::vector<size_t>> idx_data;
};



