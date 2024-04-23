from meshpy.tet import MeshInfo, build,Options
import numpy as np

def createTetraFace():
  TetraFace=np.zeros(shape=(4,3),dtype=np.int64)
  TetraFace[0,:]=[3,2,1]
  TetraFace[1,:]=[3,0,2]
  TetraFace[2,:]=[1,0,3]
  TetraFace[3,:]=[1,2,0]
  return TetraFace

def ansysMeshGen(mesh,mshF):
  TetraFace=createTetraFace()
  print(TetraFace)
  points=np.array(mesh.points)
  elements=np.array(mesh.elements)+1
  neighbors=np.array(mesh.neighbors)
  faces=np.array(mesh.faces)
  face_markers=np.array(mesh.face_markers)
  print('face_markers=',face_markers)


  TetraFace=createTetraFace()


  wall_type=[]
  interior_type=[]
  pressureOutlet_type=[]
  velocityInlet_type=[]
  nb=-1
  for j in range(0,len(neighbors)):
    # print(j,elements[j,:])
    # quit()
    for i in range(0,len(neighbors[j])):
      # if (face_markers[j]!=0):
      #   if (face_markers[j]==5):
      #     pressureOutlet_type
      
      if (neighbors[j,i]==-1):
        nb=nb+1
        i0=TetraFace[i,0]
        i1=TetraFace[i,1]
        i2=TetraFace[i,2]
        word=(str(format(elements[j,i0], 'x'))+' '+
          str(format(elements[j,i1], 'x'))+' '+
          str(format(elements[j,i2], 'x'))+' '+
          str(format(j+1, 'x'))+' '+'0\n')
        if (face_markers[nb]==3):
          wall_type.append(word)
        elif (face_markers[nb]==5):
          pressureOutlet_type.append(word)
        elif (face_markers[nb]==10):
          velocityInlet_type.append(word)
        print(nb,elements[j,i0],elements[j,i1],elements[j,i2])
      else:
        if (neighbors[j,i]>j):
          nb=nb+1
          i0=TetraFace[i,0]
          i1=TetraFace[i,1]
          i2=TetraFace[i,2]
          interior_type.append(
          str(format(elements[j,i0], 'x'))+' '+
          str(format(elements[j,i1], 'x'))+' '+
          str(format(elements[j,i2], 'x'))+' '+
          str(format(j+1, 'x'))+' '+
          str(format(neighbors[j,i]+1, 'x')+'\n'))

          print(nb,elements[j,i0],elements[j,i1],elements[j,i2])



  # print(np.array(mesh.elements))
  print(np.array(mesh.neighbors))
  # print(len(np.array(mesh.neighbors)))
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

    # wall_type,interior_type
    startnb=1
    endnb=len(interior_type)
    f.write('(13 (1 '+format(startnb, 'x')+' '+format(endnb, 'x')+' 2 3)(\n')
    for i in range(0,len(interior_type)):
      f.write(interior_type[i])
    f.write('\n')
    f.write('))\n')
    startnb=endnb+1
    endnb=endnb+len(wall_type)
    f.write('(13 (2 '+format(startnb, 'x')+' '+format(endnb, 'x')+' 3 3)(\n')
    for i in range(0,len(wall_type)):
      f.write(wall_type[i])
    f.write('))\n')
    face_gp=2
    if (len(pressureOutlet_type)>0):
      face_gp=face_gp+1
      startnb=endnb+1
      endnb=endnb+len(pressureOutlet_type)
      f.write('(13 ('+str(face_gp)+' '+format(startnb, 'x')+' '+format(endnb, 'x')+' 5 3)(\n')
      for i in range(0,len(pressureOutlet_type)):
        f.write(pressureOutlet_type[i])
      f.write('))\n')
    if (len(velocityInlet_type)>0):
      face_gp=face_gp+1
      startnb=endnb+1
      endnb=endnb+len(velocityInlet_type)
      f.write('(13 ('+str(face_gp)+' '+format(startnb, 'x')+' '+format(endnb, 'x')+' a 3)(\n')
      for i in range(0,len(velocityInlet_type)):
        f.write(velocityInlet_type[i])
      f.write('))\n')

    f.write('\n')
    f.write('))\n')

    f.write('(12 (3 1 '+format(len(elements), 'x')+' 1 2))\n')
    f.write('(45 (1 interior interior-fluid)())\n')
    f.write('(45 (2 wall wall-fluid)())\n')
    face_gp=2
    if (len(pressureOutlet_type)>0):
      face_gp=face_gp+1
      f.write('(45 ('+str(face_gp)+' pressure-outlet outlet)())\n')
    if (len(velocityInlet_type)>0):
      face_gp=face_gp+1
      f.write('(45 ('+str(face_gp)+' velocity-inlet inlet)())\n')
    face_gp=face_gp+1
    f.write('(45 ('+str(face_gp)+' fluid fluid)())\n')

#(45 (1 interior interior-fluid)())
# (45 (2 fluid fluid)())
# (45 (5 wall wall-fluid)())
# (45 (6 velocity-inlet inlet)())
# (45 (7 pressure-outlet outlet)())
# (45 (8 wall wall)())


    f.close()

'''  
f.write((0 grid written by ANSYS Meshing\n')
f.write(   nodes:       (10 (id start end type) (x y z ...))\n')
f.write(              faces:       (13 (id start end type elemType)\n')
f.write(                (v-0 v-1 .. v-n right-cell left-cell ...))\n')
f.write(   cells:       (12 (id start end type elemtype))\n')
f.write(   parent-face: (59 (start end parent child) (nchilds child0 child1 ...))\n')
f.write()\n')
'''