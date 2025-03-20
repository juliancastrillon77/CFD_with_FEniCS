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

    plt.figure(Figure, figsize=(14, 7), dpi=100)
    plt.grid(color=(147/255, 147/255, 147/255), linestyle='--', linewidth=0.2)
    plt.title(Title, color='white', fontsize=17)
    plt.ylabel(YLabel, color='white', fontsize=12)
    plt.xlabel('x', color='white', fontsize=12)
    ax = plt.gca()
    Col = (47/255, 47/255, 47/255)
    ax.set_facecolor(Col)
    plt.gcf().patch.set_facecolor(Col) 
    ax.spines['bottom'].set_color('black')
    ax.spines['left'].set_color('black')
    ax.spines['top'].set_color('black')
    ax.spines['right'].set_color('black')
    ax.tick_params(colors='white')
    plt.plot(x, Xval, color='white')

    plt.show()







