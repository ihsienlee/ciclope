from meshpy.tet import MeshInfo, build,Options
import numpy as np
import copy

def createTetraFace():
  TetraFace=np.zeros(shape=(4,3),dtype=np.int64)
  TetraFace[0,:]=[3,2,1]
  TetraFace[1,:]=[3,0,2]
  TetraFace[2,:]=[1,0,3]
  TetraFace[3,:]=[1,2,0]
  return TetraFace

# def createAnsysMsh(mesh,mshF,boundFaceInd,boundFaceType):
def createAnsysMsh(mesh,mshF):
  points=np.array(mesh.points)
  elements=np.array(mesh.elements) 
  element_attributes=np.array(mesh.element_attributes)
  faces=np.array(mesh.faces)
  face_markers=np.array(mesh.face_markers)
  mtr=np.array(mesh.point_metric_tensors)
  neighbors=np.array(mesh.neighbors)
  uniqueFaces,uniqueFaces_counts=np.unique(face_markers, return_counts=True)

  faceInd=-1
  cr=np.full(shape=len(faces),fill_value=-1,dtype=int)
  cl=np.full(shape=len(faces),fill_value=-1,dtype=int)
  for i in range(0,len(neighbors)):
    for j in range(0,4):
      if neighbors[i,j]==-1 or neighbors[i,j]>i:      
        faceInd=faceInd+1
        cr[faceInd]=i
        cl[faceInd]=neighbors[i,j]      
        # print(elements[i][TetraFace[j]]==faces[faceInd])

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
      elif  uniqueFaces[i] == 4:
        f.write('(13 (4 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' a 3)(\n')
      elif  uniqueFaces[i] == 5:
        f.write('(13 (5 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' a 3)(\n')
      elif  uniqueFaces[i] == 6:
        f.write('(13 (6 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 5 3)(\n')
      elif  uniqueFaces[i] == 7:
        f.write('(13 (7 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 3 3)(\n')
      elif  uniqueFaces[i] == 8:
        f.write('(13 (8 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 5 3)(\n')
      elif  uniqueFaces[i] == 9:
        f.write('(13 (9 '+format(startnb+1, 'x')+' '+format(startnb+len(ind), 'x')+' 3 3)(\n')
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
      elif  uniqueFaces[i] == 4:
        f.write('(45 (4 velocity-inlet x-min)())\n')
      elif  uniqueFaces[i] == 5:
        f.write('(45 (5 velocity-inlet x-max)())\n')
      elif  uniqueFaces[i] == 6:
        f.write('(45 (6 pressure-outlet y-min)())\n')
      elif  uniqueFaces[i] == 7:
        f.write('(45 (7 wall y-max)())\n')
      elif  uniqueFaces[i] == 8:
        f.write('(45 (8 pressure-outlet z-min)())\n')
      elif  uniqueFaces[i] == 9:
        f.write('(45 (9 wall z-max)())\n')
    f.close()



def ansysMeshGen(mesh,mshF):
  TetraFace=createTetraFace()
  print(TetraFace)
  points=np.array(mesh.points)
  elements=np.array(mesh.elements)+1
  neighbors=np.array(mesh.neighbors)
  faces=np.array(mesh.faces)
  face_markers=np.array(mesh.face_markers)
  # print('face_markers=',face_markers)
  # quit()

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
        # print(nb,elements[j,i0],elements[j,i1],elements[j,i2])
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

          # print(nb,elements[j,i0],elements[j,i1],elements[j,i2])



  # print(np.array(mesh.elements))
  # print(np.array(mesh.neighbors))
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