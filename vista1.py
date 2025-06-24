#%% Importar librerías
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
from matplotlib.tri import Triangulation
from matplotlib.ticker import FuncFormatter

dicmap = {  
    'Empty Channel': {
        '0.01': 'data/zero/len_15/vel_0.01/size_0.25/conc_aprox.vtk',
    },
    'Three Spacers': {
        '0.01': 'data/three/len_15/vel_0.01/size_0.25/conc_aprox.vtk/conc_aprox.vtk',
    }
}

pathToSave = "image"

def plot_concentration(case, velocity, dicmap=dicmap, pathToSave=pathToSave):
    path = dicmap[case][velocity]
    mesh = pv.read(path)

    # Extraer los puntos y las celdas (triangulos y otros)
    points = mesh.points
    triangles = mesh.cells_dict[np.uint8(5)]
    cell_data = mesh.cell_data['phih']  # Sustituir 'phih' si el nombre es diferente

    # Inicializar arrays para acumular datos y conteos
    point_data = np.zeros(len(points))    # Concentración acumulada por punto
    point_counts = np.zeros(len(points))  # Numero de contribuciones por punto

    # Asignar los valores de las celdas a los puntos (Interpolacion a Lagrange)
    for i, tri in enumerate(triangles):
        for point in tri:
            point_data[point] += cell_data[i]
            point_counts[point] += 1

    # Avoid divisions by zero
    point_counts[point_counts == 0] = 1
    point_data /= point_counts

    # Crear una triangulacion para matplotlib
    triangulation = Triangulation(points[:, 0]*1000, points[:, 1]*1000, triangles)

    #% Plot Full
    ##Filter data near 600 (makeup)
    point_data[(point_data < 600) & (point_data > 599)] = 600
    plt.figure(figsize=(20, 5), dpi=200)
    # a = [0, 0, 0]
    # b = [mesh.bounds[1], 0, 0]  # Línea horizontal en X
    # n_points = 200
    # line_data = mesh.sample_over_line(a, b, n_points)
    # vmax = line_data['phih'].max()
    levels = np.linspace(600, 800, 11)  # Niveles entre 600 y 800
    plt.tricontourf(triangulation, point_data, cmap="Spectral_r", levels=levels)
    cbar = plt.colorbar(
        orientation='vertical',
        aspect=10,
        anchor=(-0.2, 0.5),
    )
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
    cbar.set_label('mol/m3', size=15)
    cbar.ax.tick_params(labelsize=15)
    plt.xlim(0, mesh.bounds[1]*1000)
    plt.ylim(0, mesh.bounds[3]*1000)
    plt.xticks([0, 3.75, 7.5, 11.25, 15])
    plt.yticks(np.round(np.linspace(0, 0.72, 5), 2))
    plt.tick_params(labelsize=15)
    plt.ylabel("Height (mm)", fontsize=15)
    plt.xlabel("Length (mm)", fontsize=15)
    plt.title(r"Concentration Field - %s" % case, fontsize=20)
    ## Save figure
    print(f"{pathToSave}/conc_{case.replace(' ','')}_{velocity}.png")
    plt.savefig(f"{pathToSave}/conc_{case.replace(' ','')}_{velocity}.png", bbox_inches='tight')
    plt.show()
    plt.close()

    #% Plot Zoom
    range = [0.60, 0.72]
    plt.figure(figsize=(20, 5), dpi=200)
    levels = np.linspace(600, 800, 11)  # Levels between 600 y 800
    plt.tricontourf(triangulation, point_data, cmap="Spectral_r", levels=levels)
    cbar = plt.colorbar(
        orientation='vertical',
        aspect=10,
        anchor=(-0.2, 0.5),
    )
    cbar.ax.yaxis.set_major_formatter(FuncFormatter(lambda x, _: f'{x:.0f}'))
    cbar.set_label('mol/m3', size=15)
    cbar.ax.tick_params(labelsize=15)

    plt.ylim(range)
    plt.xticks([0, 3.75, 7.5, 11.25, 15])
    plt.yticks(np.round(np.linspace(range[0], range[1], 5), 2))
    plt.tick_params(labelsize=15)
    plt.ylabel("Height (mm)", fontsize=15)
    plt.xlabel("Length (mm)", fontsize=15)
    plt.title(r"Concentration Field Near the Top Membrane - %s" % case, fontsize=20)
    ## Save figure
    print(f"{pathToSave}/conc-zoom_{case.replace(' ','')}_{velocity}.png")
    plt.savefig(f"{pathToSave}/conc-zoom_{case.replace(' ','')}_{velocity}.png", bbox_inches='tight')
    plt.show()
    plt.close()

#%% Graficar
plot_concentration('Empty Channel', '0.01')
plot_concentration('Three Spacers', '0.01')
