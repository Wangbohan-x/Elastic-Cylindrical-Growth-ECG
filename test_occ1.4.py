from OCC.Core.gp import gp_Pnt, gp_Vec, gp_Ax2, gp_Dir, gp_Circ
from OCC.Core.TColgp import TColgp_HArray1OfPnt, TColgp_HArray1OfVec
from OCC.Core.TColStd import TColStd_HArray1OfBoolean
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire
from OCC.Core.BRepOffsetAPI import BRepOffsetAPI_MakePipe
from OCC.Core.Quantity import Quantity_Color, Quantity_TOC_RGB
from OCC.Display.SimpleGui import init_display
from OCC.Core.BRepMesh import BRepMesh_IncrementalMesh
from OCC.Extend.DataExchange import write_obj_file
from OCC.Core.StlAPI import StlAPI_Writer
from OCC.Core.GeomAPI import GeomAPI_Interpolate, GeomAPI_ProjectPointOnCurve
from OCC.Core.TColgp import TColgp_Array1OfPnt, TColgp_Array1OfVec
from OCC.Core.TColStd import TColStd_Array1OfReal
from OCC.Core.GeomAPI import GeomAPI_PointsToBSpline
from OCC.Core.TopoDS import TopoDS_Edge, TopoDS_Shape, TopoDS_Wire
from OCC.Core.BRep import BRep_Tool

from OCC.Core.TColStd import TColStd_HArray1OfReal
from math import sqrt
import numpy as np
import os

from OCC.Core.gp import gp_Pnt2d
from OCC.Core.Geom2dAPI import Geom2dAPI_Interpolate, Geom2dAPI_PointsToBSpline
from OCC.Core.TColgp import TColgp_HArray1OfPnt2d, TColgp_Array1OfPnt2d

from OCC.Display.SimpleGui import init_display


from OCC.Core.Geom import Geom_BSplineCurve
#display, start_display, add_menu, add_function_to_menu = init_display()

from OCC.Core.gp import gp_Pnt
from OCC.Core.BRepBuilderAPI import BRepBuilderAPI_MakeEdge, BRepBuilderAPI_MakeWire

def build_polyline(points):
    """
    根据输入点构建 OCC 折线（Wire）
    参数:
        points: List[gp_Pnt] - 输入点列表
    返回:
        TopoDS_Wire - 构建的折线
    """
    if len(points) < 2:
        raise ValueError("需要至少两个点才能生成折线。")

    wire_builder = BRepBuilderAPI_MakeWire()

    for i in range(len(points) - 1):
        edge = BRepBuilderAPI_MakeEdge(points[i], points[i + 1]).Edge()
        wire_builder.Add(edge)

    return wire_builder.Wire()

#test  scipy或者其他库============================
def normalize_points(points):
    coords = np.array([[p.X(), p.Y(), p.Z()] for p in points])
    mean = coords.mean(axis=0)
    std = coords.std(axis=0)
    norm_coords = (coords - mean) / std
    norm_points = [gp_Pnt(*c) for c in norm_coords]
    return norm_points, mean, std

def denormalize_points(norm_points, mean, std):
    # mean, std 是长度为3的 numpy 数组
    denorm_coords = [(p.X() * std[0] + mean[0],
                      p.Y() * std[1] + mean[1],
                      p.Z() * std[2] + mean[2]) for p in norm_points]
    denorm_points = [gp_Pnt(*c) for c in denorm_coords]
    return denorm_points

def build_spline_from_denormalized_points(norm_points, mean, std):
    denorm_points = denormalize_points(norm_points, mean, std)
    
    array = TColgp_Array1OfPnt(1, len(denorm_points))
    for i, pnt in enumerate(denorm_points):
        array.SetValue(i + 1, pnt)

    spline = GeomAPI_PointsToBSpline(array).Curve()
    return spline

def build_strict_3d_interpolated_spline(points):
    coords = np.array([[p.X(), p.Y(), p.Z()] for p in points])
    mean = coords.mean(axis=0)
    std = coords.std(axis=0)
    norm_coords = (coords - mean) / std
    norm_points = [gp_Pnt(*c) for c in norm_coords]

    #构建 OCC 格式的点和参数数组
    h_points = TColgp_HArray1OfPnt(1, len(norm_points))
    for i, p in enumerate(norm_points):
        h_points.SetValue(i + 1, p)

    h_params = TColStd_HArray1OfReal(1, len(norm_points))
    for i in range(len(norm_points)):
        h_params.SetValue(i + 1, i)

    #构建样条
    interpolator = GeomAPI_Interpolate(h_points, False, 1e-6)
    interpolator.Perform()
    spline_norm = interpolator.Curve()

    #提取参数
    poles = spline_norm.Poles()
    knots = spline_norm.Knots()
    multiplicities = spline_norm.Multiplicities()
    degree = spline_norm.Degree()
    is_periodic = spline_norm.IsPeriodic()

    #反归一化控制点
    new_poles = TColgp_Array1OfPnt(1, poles.Length())
    for i in range(1, poles.Length() + 1):
        p = poles.Value(i)
        x = p.X() * std[0] + mean[0]
        y = p.Y() * std[1] + mean[1]
        z = p.Z() * std[2] + mean[2]
        new_poles.SetValue(i, gp_Pnt(x, y, z))

    #构建反归一化样条
    if spline_norm.IsRational():
        weights = spline_norm.Weights()
        spline_real = Geom_BSplineCurve(
            new_poles,
            weights,
            knots,
            multiplicities,
            degree,
            is_periodic
        )
    else:
        spline_real = Geom_BSplineCurve(
            new_poles,
            knots,
            multiplicities,
            degree,
            is_periodic
        )

    return spline_real

def bspline():
    # the first bspline
    array = TColgp_Array1OfPnt2d(1, 5)
    array.SetValue(1, gp_Pnt2d(0, 0))
    array.SetValue(2, gp_Pnt2d(1, 2))
    array.SetValue(3, gp_Pnt2d(2, 3))
    array.SetValue(4, gp_Pnt2d(4, 3))
    array.SetValue(5, gp_Pnt2d(5, 5))
    bspline_1 = Geom2dAPI_PointsToBSpline(array).Curve()

    # the second one
    harray = TColgp_HArray1OfPnt2d(1, 5)
    harray.SetValue(1, gp_Pnt2d(0, 0))
    harray.SetValue(2, gp_Pnt2d(1, 2))
    harray.SetValue(3, gp_Pnt2d(2, 3))
    harray.SetValue(4, gp_Pnt2d(4, 3))
    harray.SetValue(5, gp_Pnt2d(5, 5))

    anInterpolation = Geom2dAPI_Interpolate(harray, False, 0.01)
    anInterpolation.Perform()
    bspline_2 = anInterpolation.Curve()

    harray2 = TColgp_HArray1OfPnt2d(1, 5)
    harray2.SetValue(1, gp_Pnt2d(11, 0))
    harray2.SetValue(2, gp_Pnt2d(12, 2))
    harray2.SetValue(3, gp_Pnt2d(13, 3))
    harray2.SetValue(4, gp_Pnt2d(15, 3))
    harray2.SetValue(5, gp_Pnt2d(16, 5))

    anInterpolation2 = Geom2dAPI_Interpolate(harray, True, 0.01)
    anInterpolation2.Perform()
    bspline_3 = anInterpolation2.Curve()

    for j in range(array.Lower(), array.Upper()+1):
        p = array.Value(j)
        display.DisplayShape(p, update=False)
    for j in range(harray.Lower(), harray.Upper()+1):
        p = harray.Value(j)
        display.DisplayShape(p, update=False)

    display.DisplayShape(bspline_1, update=False)
    display.DisplayShape(bspline_2, update=False, color='GREEN')
    display.DisplayShape(bspline_3, update=True, color='BLUE')


#Read points and tangents from file
def read_points_and_tangents(filename):
    points = []
    tangents = []
    with open(filename, 'r') as file:
        for line in file:
            values = line.strip().split()
            if len(values) != 6:
                continue
            try:
                x, y, z, tx, ty, tz = map(float, values)
                print(x, y, z, tx, ty, tz)
                pnt = gp_Pnt(x, y, z)
                vec = gp_Vec(tx, ty, tz)
                
                if vec.Magnitude() > 0:
                    vec.Normalize()
                points.append(pnt)
                tangents.append(vec)
            except ValueError:
                continue

    return points, tangents

#Create B-spline curve with tangents构造样条曲线（带切向量约束）
def create_spline_with_tangents(points, tangents):
    point_array = TColgp_HArray1OfPnt(1, len(points))
    tangent_array = TColgp_HArray1OfVec(1, len(tangents))
    tan_mask = TColStd_HArray1OfBoolean(1, len(tangents))
    for i, pnt in enumerate(points):
        point_array.SetValue(i + 1, pnt)
    for i, vec in enumerate(tangents):
        tangent_array.SetValue(i + 1, vec)
        tan_mask.SetValue(i + 1, True)
    interpolator = GeomAPI_Interpolate(point_array, False, 1e-6)
    interpolator.Load(tangent_array, tan_mask, True)
    interpolator.Perform()
    if not interpolator.IsDone():
        raise RuntimeError("Interpolation failed!")
    return interpolator.Curve()

def build_strict_interpolated_spline_from_points(points,tolerance = 1.0e-12):
    n = len(points)
    if n < 2:
        raise ValueError("至少需要两个点")

    # 计算三维弧长参数
    params = [0.0]
    for i in range(1, n):
        p1 = points[i-1]
        p2 = points[i]
        dist = sqrt((p2.X() - p1.X())**2 + (p2.Y() - p1.Y())**2 + (p2.Z() - p1.Z())**2)
        params.append(params[-1] + dist)

    array_points = TColgp_Array1OfPnt(1, n)
    for i, pt in enumerate(points):
        array_points.SetValue(i+1, pt)
    h_points = TColgp_HArray1OfPnt(array_points)

    array_params = TColStd_Array1OfReal(1, n)
    for i, u in enumerate(params):
        array_params.SetValue(i+1, u)
    h_params = TColStd_HArray1OfReal(array_params)

    interpolator = GeomAPI_Interpolate(h_points, h_params, False, tolerance)
    interpolator.Perform()
    if not interpolator.IsDone():
        raise RuntimeError("插值失败")

    return interpolator.Curve()


#Visualize spline, control points, and tangents可视化样条、控制点和切向量
def visualize_spline(display, spline, points, tangents):
    color1 = Quantity_Color(0.0, 0.0, 0.0, Quantity_TOC_RGB)
    display.DisplayShape(spline,color=color1, update=True)
    
    for pnt in points:
        display.DisplayShape(pnt, color="RED")
      
    for pnt, tangent in zip(points, tangents):
        end_point = gp_Pnt(pnt.X() + tangent.X(), pnt.Y() + tangent.Y(), pnt.Z() + tangent.Z())
        edge = BRepBuilderAPI_MakeEdge(pnt, end_point).Edge()
        display.DisplayShape(edge, color="GREEN")
        

#Create pipe shape along spline沿样条生成圆柱体
def create_pipe_shape( spline, start_point, start_direction, radius=0.016):
    spline_edge = BRepBuilderAPI_MakeEdge(spline).Edge()
    spline_wire = BRepBuilderAPI_MakeWire(spline_edge).Wire()

    circle_ax2 = gp_Ax2(start_point, gp_Dir(start_direction))
    circle = gp_Circ(circle_ax2, radius)
    circle_edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    circle_wire = BRepBuilderAPI_MakeWire(circle_edge).Wire()

    pipe = BRepOffsetAPI_MakePipe(spline_wire, circle_wire)
    return pipe.Shape()

def display_pipe_before_meshing(pipe_shape, display):
    display.DisplayShape(pipe_shape, update=True, color="BLUE")


def build_pipe_along_polyline(points: list[gp_Pnt], radius: float) -> TopoDS_Shape:
    # 构建折线路径
    wire_path = build_polyline(points)

    # 起点方向，用于构建扫掠圆（单位方向）
    direction = gp_Dir(points[1].X() - points[0].X(),
                       points[1].Y() - points[0].Y(),
                       points[1].Z() - points[0].Z())
    # 构造圆形截面
    circle_ax2 = gp_Ax2(points[0], direction)  # 起点位置与方向
    circle = gp_Circ(circle_ax2, radius)
    circle_edge = BRepBuilderAPI_MakeEdge(circle).Edge()
    circle_wire = BRepBuilderAPI_MakeWire(circle_edge).Wire()

    # 执行扫掠生成实体
    pipe = BRepOffsetAPI_MakePipe(wire_path, circle_wire)
    return pipe.Shape()

#Mesh and visualize shape 网格化模型（生成三角网）
def mesh_and_display(display, shape):
    mesh = BRepMesh_IncrementalMesh(shape, 0.01, True, 0.05, True)
    mesh.Perform()
    display.DisplayShape(shape, update=True)

#Export shape to OBJ
def export_to_obj(shape, filename):
    mesh = BRepMesh_IncrementalMesh(shape, 0.005, True, 0.02, True)
    mesh.Perform()
    write_obj_file(shape, filename)

#Export shape to STL
def export_to_stl(shape, filename):
    mesh = BRepMesh_IncrementalMesh(shape, 0.05, True, 0.03,True)   
    """ shape, 0.05, True, 0.04,True """
    mesh.Perform()
    writer = StlAPI_Writer()
    writer.SetASCIIMode(False)
    writer.Write(shape, filename)


def main():
    
    txt_path = "C:/Users/Shawn/Desktop/CylinderDetection/RC2/split/s4.txt"
    
    #txt_path = "points.txt"
    #spline_and_visualize(txt_path)

    
    points, tangents = read_points_and_tangents(txt_path)
    display, start_display, *_ = init_display()
    spline = create_spline_with_tangents(points, tangents)

  
    #spline = build_strict_3d_interpolated_spline(points)
    
  
    #visualize_spline(display, spline, points, tangents)
    pipe_shape = create_pipe_shape(spline, points[0], tangents[0], radius=0.016)
    #pipe_shape = build_pipe_along_polyline(points, radius=0.016)
    
    display_pipe_before_meshing(pipe_shape,display)

    
    #mesh_and_display(display, pipe_shape)
    #export_to_obj(pipe_shape, "C:/Users/Shawn/Desktop/CylinderDetection/RC/split/s1.obj")
    #export_to_stl(pipe_shape, "C:/Users/Shawn/Desktop/CylinderDetection/RC/split/s3.stl")
    
    
    white = Quantity_Color(1.0, 1.0, 1.0, Quantity_TOC_RGB)
    display.View.SetBackgroundColor(white)
    display.View.SetBgGradientColors(white, white, 0)
    display.View.Redraw()
    start_display()
 
   



def batch_process_and_export_stl(folder_path, begin, num_files):
    for i in range(begin, num_files + 1):
        txt_path = os.path.join(folder_path, f"s{i}.txt")
        stl_path = os.path.join(folder_path, f"s{i}.stl")

        txt_path1 = os.path.join(folder_path, f"r{i}.txt")
        with open(txt_path1, 'r') as f:
            lines = f.readlines()
            radius = float(lines[0].strip())


        print(f"Processing: {txt_path}")
        
        
            # 读取点和切向量
        points, tangents = read_points_and_tangents(txt_path)

            # 构建样条
        spline = create_spline_with_tangents(points, tangents)

            # 创建管道形状
        pipe_shape = create_pipe_shape(spline, points[0], tangents[0], radius=radius)

            # 导出为 STL
        export_to_stl(pipe_shape, stl_path)


if __name__ == '__main__':
    main()
    #bspline()
    #start_display()
    #batch_process_and_export_stl("C:/Users/Shawn/Desktop/CylinderDetection/RC3/split",1,12)
    #batch_process_and_export_stl("C:/Users/Shawn/Desktop/CylinderDetection/RC2/split",11,12,0.025)