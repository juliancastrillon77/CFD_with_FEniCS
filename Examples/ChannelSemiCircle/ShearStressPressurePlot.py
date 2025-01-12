# Julian Castrillon
# Shear stress & pressure plot
import os
import fenics as fe
import matplotlib.pyplot as plt
from   CfdFun import PAW

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle.xml')
DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle_facet_region.xml')
FS   = fe.FunctionSpace(Mesh, 'Lagrange', 1)
ID = 2

File = fe.File('Results/Pressure.xml')
Title = 'Pressure'
YLabel = 'Pressure'
Figure = 1
PAW(Mesh, DomainBoundaries, File, FS, ID, Figure, Title, YLabel)

File = fe.File('Results/ShearStress.xml')
Title = 'Shear Stress'
YLabel = 'Shear stress'
Figure = 2
PAW(Mesh, DomainBoundaries, File, FS, ID, Figure, Title, YLabel)

