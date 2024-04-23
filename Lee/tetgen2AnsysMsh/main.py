import pyvista as pv
import os
import numpy as np
from meshpy.tet import MeshInfo, build,Options
import vtk
import pyAnsys
import meshio




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

face_markers=[7,7,5,7,7,7,7,7,6,7,7,7]
# face_markers=np.zeros(len(faces),dtype=np.int64).tolist()
print(face_markers)

mtr=np.full(shape=len(points), fill_value=0.0, dtype=np.float64)
mtr[1]=0.5
mtr[8]=0.5
mtr[9]=0.5
mtr[12]=0.5

mtr[4]=0.5
mtr[13]=0.5
mtr[16]=0.5
mtr[17]=0.5

# print(meshConn.faces)
# faces = meshConn.faces.reshape(-1,4)[:,1:].tolist()
# face_markers=np.zeros(len(faces),dtype=np.int64).tolist()
# mtr=np.full(shape=meshConn.n_points, fill_value=0.0, dtype=np.float64)
# mtr=mtr*meshConn['CellQualityPoint']
tetgenF='test'
with open(tetgenF+'.mtr','w',encoding = 'utf-8') as f:
  f.write(str(len(mtr))+' 1\n')
  for i in range(0,len(mtr)):
    f.write(str(mtr[i])+'\n')
  f.close()

with open(tetgenF+'.node','w',encoding = 'utf-8') as f:
  f.write(str(len(points))+'  3  1  1\n')
  nb=-1
  for i in range(0,len(points)):
    nb=nb+1
    f.write(str(nb)+'  '+str(points[i][0])+'  '+str(points[i][1])+'  '+str(points[i][2])+'\n')
  f.close()

mesh_info = MeshInfo()
mesh_info.set_points(points)
mesh_info.set_facets(faces,face_markers)
mesh_info.load_mtr(tetgenF) 
switches='pqAmnfMT1e-16'
print('refinement switches=',switches)
mesh = build(mesh_info, options=Options(switches=switches))


points=np.array(mesh.points)
elements=np.array(mesh.elements)
element_attributes=np.array(mesh.element_attributes)
faces=np.array(mesh.faces)
face_markers=np.array(mesh.face_markers)
mtr=np.array(mesh.point_metric_tensors)
neighbors=np.array(mesh.neighbors)


cells=np.full((len(elements),5),4,dtype=np.int64) 
cells[:,1:5]=elements[:,0:4]
no_cells=len(cells)
cells=np.array(cells).ravel()                
celltypes = np.empty(no_cells, dtype=np.uint32)
celltypes[:] = vtk.VTK_TETRA
meshTETRA = pv.UnstructuredGrid(cells, celltypes, points)
meshTETRA=meshTETRA.compute_cell_quality(quality_measure='scaled_jacobian', null_value=-1.0)

meshTETRA['mtr']=mtr
meshTETRA['pointNb']=np.linspace(start=0,stop=meshTETRA.n_points,num=meshTETRA.n_points,endpoint=False,dtype=int)
meshTETRA['cellNb']=np.linspace(start=0,stop=meshTETRA.n_cells,num=meshTETRA.n_cells,endpoint=False,dtype=int)
meshTETRA.save('meshTETRA.vtu')


cells=np.full((len(faces),4),3,dtype=np.int64) 
cells[:,1:4]=faces[:,0:3]
no_cells=len(cells)
cells=np.array(cells).ravel()   
celltypes = np.empty(no_cells, dtype=np.uint32)
celltypes[:] = vtk.VTK_TRIANGLE
meshTRI = pv.UnstructuredGrid(cells, celltypes, points)
meshTRI['mtr']=mtr
meshTRI['face_markers']=face_markers
meshTRI.save('meshTRI.vtu')


TetraFace=pyAnsys.createTetraFace()
print(TetraFace)

faceInd=-1
print(face_markers[faceInd])
print(len(np.unique(face_markers)))
uniqueFaces,uniqueFaces_counts=np.unique(face_markers, return_counts=True)



# print(elements[0,:])
# print(faces[0:4,:])
# quit()

cr=np.full(shape=len(faces),fill_value=-1,dtype=int)
cl=np.full(shape=len(faces),fill_value=-1,dtype=int)
for i in range(0,len(neighbors)):
  for j in range(0,4):
    if neighbors[i,j]==-1 or neighbors[i,j]>i:      
      faceInd=faceInd+1
      cr[faceInd]=i
      cl[faceInd]=neighbors[i,j]      
      # print(elements[i][TetraFace[j]]==faces[faceInd])

mshF='tetgen.msh'
with open(mshF,'w',encoding = 'utf-8') as f:
  f.write('(0 grid written by ANSYS Meshing\n')
  f.write('   nodes:       (10 (id start end type) (x y z ...))\n')
  f.write('              faces:       (13 (id start end type elemType)\n')
  f.write('                (v-0 v-1 .. v-n right-cell left-cell ...))\n')
  f.write('   cells:       (12 (id start end type elemtype))\n')
  f.write('   parent-face: (59 (start end parent child) (nchilds child0 child1 ...))\n')
  f.write(')\n')

  f.write('(2 3)\n')
  f.write('(10 (0 1 '+format(len(points), 'x')+' 0))\n')
  f.write('(13 (0 1 '+format(len(faces), 'x')+' 0))\n')
  f.write('(12 (0 1 '+format(len(elements), 'x')+' 0))\n')

  f.write('(10 (3 1 '+format(len(points), 'x')+' 1 3)(\n')
  for i in range(0,len(points)):
    f.write(str(points[i,0])+' '+str(points[i,1])+' '+str(points[i,2])+'\n')
  f.write('))\n')

  # uniqueFaces,uniqueFaces_counts
  startnb=0
  for i in range(0,len(uniqueFaces)):
    ind=np.array(list(zip(*np.where(face_markers==uniqueFaces[i]))))[:,0]
    print("startnb,len(ind)",startnb,len(ind))
    if uniqueFaces[i] == 0:      
      # print(len(ind),uniqueFaces_counts[i])
      f.write('(13 (1 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 2 3)(\n')
    elif  uniqueFaces[i] == 5:
      f.write('(13 (5 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' a 3)(\n')
    elif  uniqueFaces[i] == 6:
      f.write('(13 (6 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 5 3)(\n')
    elif  uniqueFaces[i] == 7:
      f.write('(13 (7 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 3 3)(\n')
    for j in range(0,len(ind)):
      f.write(
        str(format(faces[ind[j]][0]+1, 'x'))+' '+
        str(format(faces[ind[j]][1]+1, 'x'))+' '+
        str(format(faces[ind[j]][2]+1, 'x'))+' '+
        str(format(cr[ind[j]]+1, 'x'))+' '+
        str(format(cl[ind[j]]+1, 'x'))+'\n'
        )
    f.write('\n')
    f.write('))\n')
    startnb=startnb+len(ind)

  f.write('(12 (2 1 '+format(len(elements), 'x')+' 1 2))\n')
  for i in range(0,len(uniqueFaces)):  
    if uniqueFaces[i] == 0: 
      f.write('(45 (1 interior interior-fluid)())\n')
      f.write('(45 (2 fluid fluid)())\n')
    elif  uniqueFaces[i] == 5:
      f.write('(45 (5 velocity-inlet inlet)())\n')
    elif  uniqueFaces[i] == 6:
      f.write('(45 (6 pressure-outlet outlet)())\n')
    elif  uniqueFaces[i] == 7:
      f.write('(45 (7 wall wall)())\n')
  f.close()

print("len(elements)", len(elements))
print("len(faces)", len(faces))
print("len(points)", len(points))

ws=os.getcwd()

filename=os.path.join(ws,'tetgen.msh')
# filename=os.path.join(ws,'../readMsh/ANSYS_MSH_example','bcexample.msh')

mesh = meshio.read(
    filename,  # string, os.PathLike, or a buffer/open file
    file_format="ansys",  # optional if filename is a path; inferred from extension
    # see meshio-convert -h for all possible formats
)
print(mesh)
# print(mesh.points)
# print(dir(mesh))
# print(mesh.cell_data_to_sets)
# print(mesh.cell_sets_to_data)
# print(mesh.cells)
# print(mesh.cells_dict['triangle'])
# print(len(mesh.cells_dict['triangle']))


