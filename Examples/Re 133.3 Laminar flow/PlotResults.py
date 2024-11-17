# Julian Castrillon
# Shear stress & pressure plot
import os
import fenics as fe
import numpy as np
import matplotlib.pyplot as plt

def clear_console():
    os.system('cls' if os.name == 'nt' else 'clear')
clear_console()

Mesh = fe.Mesh('Mesh/2D Mesh/2DChannelSemiCircle.xml')
DomainBoundaries = fe.MeshFunction('size_t', Mesh, 'Mesh/2D Mesh/2DChannelSemiCircle_facet_region.xml')

Vel  = fe.VectorElement('Lagrange', Mesh.ufl_cell(), 2)
Pres = fe.FiniteElement('Lagrange', Mesh.ufl_cell(), 1)
M    = fe.MixedElement([Vel, Pres])
FS   = fe.FunctionSpace(Mesh, M)
SSFS = fe.TensorFunctionSpace(Mesh, 'Lagrange', 1)

u = fe.Function(FS.sub(0).collapse())  # Function for velocity
p = fe.Function(FS.sub(1).collapse())  # Function for pressure
ShearStress = fe.Function(SSFS)

fe.File('Results/Velocity.xml') >> u
fe.File('Results/Pressure.xml') >> p
fe.File('Results/ShearStress.xml') >> ShearStress

Boundarypoints = []
for Facet in fe.facets(Mesh):
    if DomainBoundaries[Facet] == 2:
        for Vertex in fe.vertices(Facet):
            Point = Vertex.point()
            Boundarypoints.append((Point.x(), Point.y()))

PointsList = list(set(Boundarypoints))
x, y = zip(*PointsList)  # Unpack x and y values

PressureVec = []
ShearStressVec = []
PointsWithPressureAndShear = []

for point in PointsList:
    x, y = point
    PressureXY    = p([x, y])  # Evaluate shear at (x, y)
    ShearStressXY = ShearStress([x, y])  # Evaluate shear at (x, y)
    ShearStressXY = np.sqrt(np.sum(ShearStressXY**2))
    PressureVec.append(PressureXY)
    ShearStressVec.append(ShearStressXY)
    PointsWithPressureAndShear.append((x, y, PressureXY, ShearStressXY))

SortedPoints = sorted(PointsWithPressureAndShear, key=lambda point: point[0])
x, y_sorted, yP, yS = zip(*SortedPoints)

plt.figure(1)
plt.grid()
plt.title('Shear stress bottom wall')
plt.ylabel('y')
plt.xlabel('x')
plt.plot(x, yS, linestyle='solid', color='#4DBEEE')

plt.figure(2)
plt.grid()
plt.title('Pressure bottom wall')
plt.ylabel('y')
plt.xlabel('x')
plt.plot(x, yP, linestyle='solid', color='#4DBEEE')

plt.show()








