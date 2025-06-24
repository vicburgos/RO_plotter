#%% Importar librerías
import pyvista as pv
import matplotlib.pyplot as plt
import numpy as np
from scipy.interpolate import griddata
import numpy as np
import matplotlib.tri as tri

dicmap = {  
    'Empty Channel': {
        '0.01': '/home/vburgos/RO_operation/data/zero/len_15/vel_0.01/size_0.25/flux_aprox.vtk',
    },
    'Three Spacers': {
        '0.01': '/home/vburgos/RO_operation/data/three/len_15/vel_0.01/size_0.25/flux_aprox.vtk',
    }
}

pathToSave = "image"

def plot_velocity(case, velocity, dicmap=dicmap, pathToSave=pathToSave):
    print(f"Plotting {case} - {velocity}")
    ## Load data
    file_path = dicmap[case][velocity]
    mesh = pv.read(file_path)

    points = mesh.points
    triangles = mesh.cells_dict[np.uint8(5)]
    cell_data_u = mesh.cell_data['uh'][:,0]
    cell_data_v = mesh.cell_data['uh'][:,1]

    point_data_u = np.zeros(len(points))
    point_data_v = np.zeros(len(points))
    point_counts = np.zeros(len(points))

    for i, triangulo in enumerate(triangles):
        for point in triangulo:
            point_data_u[point] += cell_data_u[i]
            point_data_v[point] += cell_data_v[i]
            point_counts[point] += 1

    point_counts[point_counts == 0] = 1
    point_data_u /= point_counts
    point_data_v /= point_counts

    triangulation = tri.Triangulation(points[:, 0]*1000, points[:, 1]*1000, triangles)

    ## Plot velocity
    if velocity == "0.01":
        ticks_velocity = np.round(np.linspace(0, 2.77*0.01, 4), 2)
    else:
        ticks_velocity = np.round(np.linspace(0, 0.18, 4), 2)
    plt.figure(figsize=(20, 5), dpi=200)
    plt.tripcolor(
        triangulation, 
        np.sqrt(point_data_u**2+point_data_v**2), 
        cmap="viridis",
        vmin=0,
        vmax=ticks_velocity[-1],
        shading='gouraud'
    )
    cbar=plt.colorbar(
        orientation='vertical', 
        aspect=10,
        anchor=(-0.2, 0.5),
    )
    cbar.set_label('m/s', size=15)
    cbar.ax.tick_params(labelsize=15)
    cbar.set_ticks(ticks_velocity)

    # Agregar las streamlines
    grid_x, grid_y = np.meshgrid(
        np.linspace(0, points[:, 0].max(), 200),
        np.linspace(0, points[:, 1].max(), 20)
    )
    number_seed = 8
    start_x = np.ones(number_seed)*points[:, 0].max()*0.62
    start_y = np.linspace(
        points[:, 1].max() * 0.15,
        points[:, 1].max() * 0.85,
        number_seed
    )
    start_points = np.array([start_x, start_y]).T  # Crear las posiciones de las semillas

    grid_u = griddata(points[:, :2], point_data_u, (grid_x, grid_y), method='linear')
    grid_v = griddata(points[:, :2], point_data_v, (grid_x, grid_y), method='linear')

    grid_u = np.nan_to_num(grid_u, nan=0.0)
    grid_v = np.nan_to_num(grid_v, nan=0.0)
    plt.streamplot(
        grid_x * 1000, grid_y * 1000,  # Convertir las coordenadas a mm
        grid_u, grid_v,
        color='white', linewidth=2, density=0.4, arrowsize=0,
        broken_streamlines = False,
        start_points=start_points * 1000
    )
    plt.tick_params(labelsize=15)
    plt.xticks([0, 3.75, 7.50, 11.25, 15])
    plt.yticks(np.round(np.linspace(0, 0.72, 5), 2))
    plt.ylabel("Height (mm)", fontsize=15)
    plt.xlabel("Length (mm)", fontsize=15)
    #Numero de ticks en el eje y: 5
    # plt.yticks(np.round(np.linspace(0, points[:, 1].max()*1000, 5), 2))
    plt.title(r"Velocity Field - %s" % case, fontsize=20)
    ## Save figure
    plt.show()
    print(f"{pathToSave}/velocity_{case.replace(' ','')}_{velocity}.png")
    plt.savefig(f"{pathToSave}/velocity_{case.replace(' ','')}_{velocity}.png", bbox_inches='tight', dpi=200)
    plt.close()

    ## Plot zoom
    plt.figure(figsize=(20, 5), dpi=200)
    y_range = 0.94
    grid_x, grid_y = np.meshgrid(
        np.linspace(points[:, 0].min(),         points[:, 0].max(), 200),
        np.linspace(points[:, 1].max()*y_range, points[:, 1].max(), 20)
    )

    grid_u = griddata(points[:, :2], point_data_u, (grid_x, grid_y), method='linear')
    grid_v = griddata(points[:, :2], point_data_v, (grid_x, grid_y), method='linear')

    grid_u = np.nan_to_num(grid_u, nan=0.0)
    grid_v = np.nan_to_num(grid_v, nan=0.0)

    plt.tripcolor(
        triangulation, 
        np.sqrt(point_data_u**2+point_data_v**2), 
        vmin=0, 
        vmax=ticks_velocity[-1],
        shading='gouraud'
    )
    cbar=plt.colorbar(
        orientation='vertical', 
        aspect=10,
        anchor=(-0.2, 0.5),
    )
    cbar.set_label('m/s', size=15)
    cbar.ax.tick_params(labelsize=15)
    cbar.set_ticks(ticks_velocity)

    plt.streamplot(
        grid_x * 1000, grid_y * 1000,  # Convertir las coordenadas a mm
        grid_u, grid_v,
        color='white', linewidth=2, density=1.3, arrowsize=1,
        broken_streamlines = True
    )
    plt.xlim(0, points[:, 0].max()*1000)
    plt.ylim(0.69, points[:, 1].max()*1000)
    plt.xticks([0, 3.75, 7.50, 11.25, 15])
    plt.yticks(np.round(np.linspace(0.69, 0.72, 4), 2))
    plt.tick_params(labelsize=15)
    plt.ylabel("Height (mm)", fontsize=15)
    plt.xlabel("Length (mm)", fontsize=15)
    plt.title(r"Velocity Field Near the Top Membrane - %s" % case, fontsize=20)
    ## Save figure
    print(f"{pathToSave}/velocity-zoom_{case.replace(' ','')}_{velocity}.png")
    plt.savefig(f"{pathToSave}/velocity-zoom_{case.replace(' ','')}_{velocity}.png", bbox_inches='tight', dpi=200)
    plt.show()
    plt.close()


#%% Graficar
plot_velocity('Empty Channel', '0.01')
plot_velocity('Three Spacers', '0.01')