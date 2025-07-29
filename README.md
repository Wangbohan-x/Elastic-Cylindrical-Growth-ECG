# Elastic-Cylindrical-Growth-ECG
主要依赖库：
c++部分（点云读取，钢筋提取，参数输出）依赖库均用vcpkg安装
cgal库-圆柱体拟合ransac算法
CGAL 5.5.2
依赖（vcpkg 会自动安装）：
    Boost
    GMP, MPFR
    Eigen3
具体安装步骤参照官网：
https://doc.cgal.org/latest/Manual/windows.html#title0
	Using CGAL on Windows (with Visual C++) 
		Installing CGAL with the Vcpkg Library Manager

ceres库-圆柱体最小二乘法参数优化
依赖（vcpkg 会自动安装）：
    Eigen3
    glog
    CXSparse / SuiteSparse（可选）
具体安装步骤参照官网：
http://ceres-solver.org/installation.html#using-a-library-manager
	Windows
		Using a Library Manager

python部分（建模，输出stl）
occ库-样条拟合，模型输出
pythonocc-core=7.6.1（对应 OCCT 7.6.0）
使用 conda安装
参照论坛内容：
https://blog.csdn.net/gitblog_09321/article/details/142223861
