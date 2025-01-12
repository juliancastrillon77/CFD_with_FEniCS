# Julian Castrillon
# CFD Functions

import fenics as fe
import matplotlib.pyplot as plt

def PAW(Mesh, DomainBoundaries, File, FS, ID, Figure, Title, YLabel): # Plot along wall

    X = fe.Function(FS)
    File >> X

    Boundarypoints = []
    for Facet in fe.facets(Mesh):
        if DomainBoundaries[Facet] == ID:
            for Vertex in fe.vertices(Facet):
                Point = Vertex.point()
                S = Point.x()
                Boundarypoints.append((Point.x(), Point.y()))

    PointsSort = sorted(Boundarypoints, key=lambda point: point[0])

    x, y = zip(*PointsSort)
    Xval = [X(xi, yi) for xi, yi in zip(x, y)]  # Evaluating X at each point

    plt.figure(Figure, figsize=(15, 7), dpi=100)  # Example size
    plt.grid()
    plt.title(Title)
    plt.ylabel(YLabel)
    plt.xlabel('x')
    plt.plot(x, Xval, color='b')
    plt.show()







