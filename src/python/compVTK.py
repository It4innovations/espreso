#####!/opt/paraview/ParaView-5.0.1-Qt4-OpenGL2-MPI-Linux-64bit/bin/pvpython
# CHANGE FIRST LINE ACCORDING TO THE LOCATION OF YOUR PVPYTHON ON YOUR MACHINE

typeToSave = 'onlySurf'  # 'full', 'onlySurf' , 'onlySurfAndReduc3angles'

from mpi4py import MPI
from paraview.simple import *
import glob
#import os.path
import os
import errno
import math
#
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#
comm = MPI.COMM_WORLD
rank = comm.Get_rank()
size = comm.Get_size() 
#
path = 'OUT'
vtk_files = 'result'
# SET SUFFIX (VTU: WHOLE VOLUME, VTP: SURFACE ONLY)
if typeToSave == 'full':
  suff = 'vtu'
elif typeToSave == 'onlySurf' or typeToSave=='onlySurfAndReduc3angles':
  suff = 'vtp'
# EXISTING DIRECTORY CHECKING
try:
  os.makedirs(path)
  os.system('OUT/*')
except OSError as exception:
  if exception.errno != errno.EEXIST:
    raise
# FILES SAVING FUNCTION
def saveFile(i_file,rank):
  r = OpenDataFile(i_file) 
  _in = int(filter(str.isdigit,i_file))
  if typeToSave=='full':
    tmp = r
  elif typeToSave == 'onlySurf':
    dr = ExtractSurface(Input=r)
    tmp = Decimate(Input=dr)
    #tmp = dr
  elif typeToSave == 'onlySurfAndReduc3angles':
    dr= Delaunay3D(Input=r)
    s = ExtractSurface(Input=dr)
    tmp = Decimate(Input=s)
  SaveData(path+'/output'+str(_in)+'.'+suff, proxy=tmp)
  return 0
#
#------------create vtm file --------------------------------------------------
# LIST_: LIST OF ALL VTK FILES IN CURRENT DIRECTORY GOING TO BE CONVERTED

comm.Barrier()

list_ = glob.glob(vtk_files+'*.vtk')
_n_list = len(list_)
# _DEL : NUMB. OF VTK FILES, EACH MPI IS PROCESSING
_del = int(math.ceil(float(_n_list) / size))
if rank==0:
  print '_del:  ', _del
  print 'size:  ', size
  print '_n_list:',_n_list 
#
for i in range(_n_list):
  if i >= rank*_del and i < (rank+1)*_del:
    _in = int(filter(str.isdigit,list_[i]))
    print 'rank',rank,'_in',_in
    saveFile(list_[i],rank)
#
if size * _del < _n_list  and rank == 0:
  for i in range(size*_del,_n_list):
    _in = int(filter(str.isdigit,list_[i]))
    saveFile(list_[i],rank)
    print 'rank',rank,'_in',_in, '++++++++++++++++++'
#

comm.Barrier()

if rank == 0:
  list_ = glob.glob(path+'/*.'+suff)
  print list_
#
  f = open('set.vtm', 'w')
  str0 = '<VTKFile type="vtkMultiBlockDataSet" version="1.0" byte_order="LittleEndian" header_type="UInt64" compressor="vtkZLibDataCompressor">\n'
  f.write(str0)
  str1 = '  <vtkMultiBlockDataSet>\n'
  f.write(str1)
#
  for i in range(len(list_)):
    i_file = list_[i]
    f.write('    <DataSet index="'+str(i)+'" file="'+i_file+'">\n')
    f.write('    </DataSet>\n')
#
  str2 = '  </vtkMultiBlockDataSet>\n'
  f.write(str2)
  str3 = '</VTKFile>\n'
  f.write(str3)
#
  f.close()

