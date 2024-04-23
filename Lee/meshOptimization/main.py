import pyvista as pv
import os
import numpy as np
from meshpy.tet import MeshInfo, build,Options
import vtk
import pyAnsys

pv.global_theme.transparent_background = True

path=os.getcwd()

def findBoundInd(mesh,bounds):
  meshConn_cell_centers=meshConn.cell_centers()
  meshConn_cell_centers.save("meshConn_cell_centers.vtk")

  minx=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[0])
  dis=abs(minx-meshConn_cell_centers.points[:,0])
  minXind=dis<1.e-2
  maxx=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[1])
  dis=abs(maxx-meshConn_cell_centers.points[:,0])
  maxXind=dis<1.e-2
  
  miny=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[2])
  dis=abs(miny-meshConn_cell_centers.points[:,1])
  minYind=dis<1.e-2
  maxy=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[3])
  dis=abs(maxy-meshConn_cell_centers.points[:,1])
  maxYind=dis<1.e-2
  
  minz=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[4])
  dis=abs(minz-meshConn_cell_centers.points[:,2])
  minZind=dis<1.e-2

  maxz=np.full(shape=meshConn_cell_centers.n_points,fill_value=bounds[5])
  dis=abs(maxz-meshConn_cell_centers.points[:,2])
  maxZind=dis<1.e-2

  return minXind,maxXind,minYind,maxYind,minZind,maxZind

case=0
if case==1:
  dirPath = os.path.join(path,'structure')
  results = [f for f in os.listdir(dirPath) if os.path.isfile(os.path.join(dirPath, f))]
  outputdirPath = os.path.join(path,'result')
  for i in range(0,len(results)):
    filename=os.path.join(dirPath,results[i])
    mesh=pv.read(filename)
    meshConn=mesh.connectivity(extraction_mode='largest', variable_input=None, scalar_range=None, scalars=None, label_regions=True, region_ids=None, point_ids=None, cell_ids=None, closest_point=None,)
    # meshConn=meshConn.compute_cell_quality(quality_measure='scaled_jacobian', null_value=-1.0)
    output_filename=os.path.join(outputdirPath,results[i])
    meshConn.save(output_filename)

min_mtr=0.005
# min_mtr=0.011709942281448605
# min_mtr=0.0001562256705241939
case=0
if case==0:
  # sourceF='0_51'
  # sourceF='50_51'
  sourceF='100_51'
  sourceF='150_51'
  sourceF='200_51'
  sourceF='250_51'
  sourceF='300_51'
  sourceF='350_51'
  sourceF='400_51'
  sourceF='450_51'
  sourceF='499_51'
  filename=os.path.join(path,'structure','result_modified',sourceF+'.stl')
  mesh=pv.read(filename)



  meshConn=mesh.connectivity(extraction_mode='largest', variable_input=None, scalar_range=None, scalars=None, label_regions=True, region_ids=None, point_ids=None, cell_ids=None, closest_point=None,)

  bounds=meshConn.bounds
  # print(meshConn)
  print("sourceF bounds",bounds)
  # print(type(bounds[0]))
  # quit()
  minXind,maxXind,minYind,maxYind,minZind,maxZind=findBoundInd(meshConn,bounds) 

  # meshConn.save('0_51_Conn.stl')

  smooth = meshConn.smooth()
  # meshConn.save('0_51_ConnSmooth.stl')
  smooth1000 = meshConn.smooth(n_iter=1000)
  # meshConn.save('0_51_ConnSmooth1000.stl')

  smooth_w_taubin = meshConn.smooth_taubin(n_iter=2, pass_band=0.1)
  # smooth_w_taubin.save('0_51_ConnTaubin.stl')
  smooth_w_taubin50 = meshConn.smooth_taubin(n_iter=50, pass_band=0.05)
  # smooth_w_taubin.save('0_51_ConnTaubin50.stl')


  # meshConn=smooth_w_taubin

  meshConn=meshConn.compute_cell_quality(quality_measure='min_angle', null_value=-1.0)
  # print(np.nanmax(meshConn['CellQuality']))
  # print(np.nanmin(meshConn['CellQuality']))
  # # meshConn.save('meshConn.vtk')
  # quit()



  print("meshConn.n_points=",meshConn.n_points)
  print("max faces=",np.max(meshConn.faces))

  faces = meshConn.faces.reshape(-1,4)[:,1:]
  face_markers=np.zeros(len(faces),dtype=np.int64)
  boundFaceInd=np.array([4,5,6,7,8,9])
  boundFaceType=np.array([4,5,7,7,7,7])
  face_markers[minXind]=boundFaceInd[0]
  face_markers[maxXind]=boundFaceInd[1]
  face_markers[minYind]=boundFaceInd[2]
  face_markers[maxYind]=boundFaceInd[3]
  face_markers[minZind]=boundFaceInd[4]
  face_markers[maxZind]=boundFaceInd[5]
  meshConn['face_markers']=face_markers
  face_markers=face_markers.tolist()

  meshConn=meshConn.compute_cell_sizes()
  print(meshConn.array_names)
  # print(meshConn['Length'])
  mtr=np.full(shape=meshConn.n_points, fill_value=3000.0, dtype=np.float64)
  meshConn['CellQualityPoint']=np.full(shape=meshConn.n_points,fill_value=90.0,dtype=np.float64)

  nb=0
  CellQualityValue=1.5
  case=0
  if case==1:
    angle2refineCheck=True
    if angle2refineCheck:
      for i in range(0,meshConn.n_cells):
        for point_id in faces[i]:        
          if meshConn['CellQuality'][i]>60.0:
            angle=12.0-meshConn['CellQuality'][i]
            if meshConn['CellQualityPoint'][point_id]>angle:
              meshConn['CellQualityPoint'][point_id]=angle

          else:
            if meshConn['CellQualityPoint'][point_id]>meshConn['CellQuality'][i]:
              meshConn['CellQualityPoint'][point_id]=meshConn['CellQuality'][i]
      meshConn['CellQualityPoint']=meshConn['CellQualityPoint']/np.float64(60.0)
      mtr=mtr*meshConn['CellQualityPoint']
      # print(mtr)
      # print(np.max(mtr), np.min(mtr))
    else:
      for i in range(0,meshConn.n_cells):
        for point_id in faces[i]: 
          if meshConn['CellQuality'][i]>CellQualityValue:
            cell = meshConn.get_cell(i)
            # print("i=",i,meshConn['Length'][i])
            cell_point_ids = cell.point_ids
            # print('cell_point_ids=',cell_point_ids)
            n_edges=cell.n_edges
            disarray=np.zeros(n_edges,dtype=float)
            for j in range(0,n_edges):
              edges = cell.get_edge(j)
              edge_point_ids = edges.point_ids
              # print('edge_point_ids=',edge_point_ids)
              a=np.array(meshConn.points[edge_point_ids[0]])
              b=np.array(meshConn.points[edge_point_ids[1]])
              disarray[j]=np.linalg.norm(a-b)
            mindis=np.nanmin(disarray)
            for j in range(0,len(cell_point_ids)):
              if mtr[cell_point_ids[j]]>mindis*CellQualityValue:
                mtr[cell_point_ids[j]]=mindis
  meshConn.save('meshConn.vtk')
  # quit()

  points=tuple(map(tuple, meshConn.points))
  print(meshConn.faces)
  faces = meshConn.faces.reshape(-1,4)[:,1:].tolist()
  # face_markers=np.zeros(len(faces),dtype=np.int64).tolist()
  print('np.min(mtr)=',np.min(mtr))
  if np.nanmin(mtr)<min_mtr:
    min_mtr_ind=np.array(list(zip(*np.where(mtr<min_mtr))))[:,0]
    mtr[min_mtr_ind]=min_mtr
    # quit()

  tetgenF='test'
  with open(tetgenF+'.mtr','w',encoding = 'utf-8') as f:
    f.write(str(len(mtr))+' 1\n')
    for i in range(0,len(mtr)):
      f.write(str(mtr[i])+'\n')
    f.close()

  # with open(tetgenF+'.node','w',encoding = 'utf-8') as f:
  #   f.write(str(len(points))+'  3  1  1\n')
  #   nb=-1
  #   for i in range(0,len(points)):
  #     nb=nb+1
  #     f.write(str(nb)+'  '+str(points[i][0])+'  '+str(points[i][1])+'  '+str(points[i][2])+'\n')
  #   f.close()

elif case==1:
  points=np.array(
    [
      [0.0, 0.0, 2.0], # 0 
      [2.0, 0.0, 2.0], 
      [2.0, 2.0, 2.0], 
      [0.0, 2.0, 2.0],
      [0.0, 0.0, 0.0], 
      [2.0, 0.0, 0.0], 
      [2.0, 2.0, 0.0], 
      [0.0, 2.0, 0.0],

      [1.0, 0.0, 2.0],
      [2.0, 1.0, 2.0],
      [1.0, 2.0, 2.0],
      [0.0, 1.0, 2.0],
      [1.0, 1.0, 2.0],

      [1.0, 0.0, 0.0],
      [2.0, 1.0, 0.0],
      [1.0, 2.0, 0.0],
      [0.0, 1.0, 0.0],
      [1.0, 1.0, 0.0],
    ])


  points=tuple(map(tuple, points))

  faces=[
      [0, 8, 12, 11],
      [11, 12, 10, 3],
      [8, 1, 9, 12],
      [12, 9, 2, 10],
      [0, 4, 13, 5, 1, 8],
      [1, 5, 14, 6, 2, 9],
      [2, 6, 15, 7, 3, 10],
      [3, 7, 16, 4, 0, 11],
      [4, 13, 17, 16],
      [16, 17, 15, 7],
      [13, 5, 14, 17],
      [17, 14, 6, 15],
    ]

  face_markers=[8,8,8,8,6,5,7,4,9,9,9,9]

  mtr=np.full(shape=len(points), fill_value=0.0, dtype=np.float64)
  mtr[1]=0.5
  mtr[8]=0.5
  mtr[9]=0.5
  mtr[12]=0.5

  mtr[4]=0.5
  mtr[13]=0.5
  mtr[16]=0.5
  mtr[17]=0.5

  tetgenF='test'
  with open(tetgenF+'.mtr','w',encoding = 'utf-8') as f:
    f.write(str(len(mtr))+' 1\n')
    for i in range(0,len(mtr)):
      f.write(str(mtr[i])+'\n')
    f.close()

mesh_info = MeshInfo()
mesh_info.set_points(points)
mesh_info.set_facets(faces,face_markers)
mesh_info.load_mtr(tetgenF) 
mesh_info.save_poly(tetgenF)
mesh_info.save_nodes(tetgenF)

switches='pqAmnfMT1e-16'
print('refinement switches=',switches)
kwargs='voroout'
mesh = build(mesh_info, options=Options(switches=switches))

# mshF='meshConn.msh'
# pyAnsys.ansysMeshGen(mesh,mshF) #old version


points=np.array(mesh.points)
elements=np.array(mesh.elements) 
element_attributes=np.array(mesh.element_attributes)
faces=np.array(mesh.faces)
face_markers=np.array(mesh.face_markers)
mtr=np.array(mesh.point_metric_tensors)
voroout=np.array(mesh.voroout)
print(dir(mesh))
print(np.nanmax(mtr))
print(np.nanmin(mtr))



cells=np.full((len(elements),5),4,dtype=np.int64) 
cells[:,1:5]=elements[:,0:4]
no_cells=len(cells)
cells=np.array(cells).ravel()                
celltypes = np.empty(no_cells, dtype=np.uint32)
celltypes[:] = vtk.VTK_TETRA
meshTETRA = pv.UnstructuredGrid(cells, celltypes, points)
meshTETRA=meshTETRA.compute_cell_quality(quality_measure='min_angle', null_value=-1.0)
meshTETRA['mtr']=mtr
meshTETRA.save('meshTETRA.vtu')

cells=np.full((len(faces),4),3,dtype=np.int64) 
cells[:,1:4]=faces[:,0:3]
no_cells=len(cells)
cells=np.array(cells).ravel()                
celltypes = np.empty(no_cells, dtype=np.uint32)
celltypes[:] = vtk.VTK_TRIANGLE
meshTRI = pv.UnstructuredGrid(cells, celltypes, points)
meshTRI=meshTRI.compute_cell_quality(quality_measure='min_angle', null_value=-1.0)

meshTRI['mtr']=mtr
meshTRI['face_markers']=face_markers
meshTRI.save('meshTRI.vtu')


quit()

mshF=sourceF+'.msh'
# pyAnsys.createAnsysMsh(mesh,mshF,boundFaceInd,boundFaceType)
pyAnsys.createAnsysMsh(mesh,mshF)
# meshTETRA_surface=meshTETRA.extract_surface()
# meshTETRA_surface.save('ConnRefinment.stl')
# meshTETRA_surface=meshTETRA_surface.compute_cell_quality(quality_measure='scaled_jacobian', null_value=-1.0)
# meshTETRA_surface.save('meshTETRA_surface.vtk')
