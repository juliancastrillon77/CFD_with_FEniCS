# Julian Castrillon
# Shear stress & pressure plot
import os
import fenics as fe
import matplotlib.pyplot as plt
from   CfdFun import PAW

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle2.xml')
DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle2_facet_region.xml')
FS = fe.FunctionSpace(Mesh, 'Lagrange', 1)
ID = 2

File = fe.File('Results/Pressure.xml') # 25087
Title = 'Pressure along bottom wall (channel semicircle example)'
YLabel = 'Gauge pressure [Pa]'
Figure = 1
PAW(Mesh, DomainBoundaries, File, FS, ID, Figure, Title, YLabel)

File = fe.File('Results/ShearStress.xml')
Title = 'Shear stress along bottom wall (channel semicircle example)'
YLabel = 'Shear stress [Pa]'
Figure = 2
PAW(Mesh, DomainBoundaries, File, FS, ID, Figure, Title, YLabel)

