## Function to read a VTK file and generate visualisations for Stream functions and Vorticity
# %% 
import numpy as np
import matplotlib.pyplot as plt
import pandas as pd
import meshio
import sys
import os
# %%
## main function
if __name__ == "__main__":

    # read the command line argument for vtk file name and exit if not provided
    
    if len(sys.argv) < 2:
        print("Please provide the vtk file name")
        sys.exit(0)
    else:
        vtk_file_name = sys.argv[1]

# %%    


    # vtk_file_name = "Mean_NRealisations_50_00000.vtk"

    ## read the vtk file using meshio    
    mesh = meshio.read(vtk_file_name, file_format="vtk")

    ## Extract the points and solution data from the mesh
    points = mesh.points
    point_data = mesh.point_data

    ## N_DOF 
    N_DOF = len(points)

    ## N
    N = int(np.sqrt(N_DOF))

    # Create a Matrix of size N_DOF x N_DOF
    U = np.zeros((N, N))
    V = np.zeros((N, N))


    ## Assign Boundary Conditions to the  Matrix U
    for i in range(N):
        xVal = i/(N-1)
        for j in range(N):
            yVal = j/(N-1)        
            if( np.isclose(xVal, 0) or np.isclose(xVal, 1) or np.isclose(yVal, 0) or np.isclose(yVal, 1)):
                U[i,j] = 0
            if (np.isclose)(yVal , 1.0):
                U[i,j] = 1.0

    ## Assign Boundary Conditions to the  Matrix V
    for i in range(N):
        xVal = i/(N-1)
        for j in range(N):
            yVal = j/(N-1)        
            if( np.isclose(xVal, 0) or np.isclose(xVal, 1) or np.isclose(yVal, 0) or np.isclose(yVal, 1)):
                V[i,j] = 0
            
    ## Create a mapping for the points coordinates from VTK to the index of the matrix A
    point_map = {}
    for i in range(N_DOF):
        point_map[round(points[i][0],6), round(points[i][1],6)] = i

    ## Assign the values of the solution to the matrix A
    for i in range(1,N-1):
        xVal = i/(N-1)
        for j in range(1,N-1):
            yVal = j/(N-1)
            U[i,j] = point_data["U_Mean1"][point_map[round(xVal,6), round(yVal,6)]] 
            V[i,j] = point_data["U_Mean2"][point_map[round(xVal,6), round(yVal,6)]] 

        
    ## Plot the Velocities U and V in contour plots with contour levels 90 and colormaps jet 
    ## plot it using ax.contourf
    fig, ax = plt.subplots(1,2, figsize=(10,5))
    ax[0].contourf(U.T, 40, cmap='jet')
    ax[1].contourf(V.T, 40, cmap='jet')
    ax[0].set_title("U")
    ax[1].set_title("V")
    plt.savefig("U_V.png",dpi=300)


    ## Compute the voriticity on the interior points
    Vorticity = np.zeros((N,N))
    for i in range(1,N-1):
        for j in range(1,N-1):
            Vorticity[i,j] =  (V[i,j+1] - V[i,j-1]) - (U[i+1,j] - U[i-1,j])

    ## Plot the Vorticity in a contour plot with contour levels 90 and colormap jet
    ## plot it using ax.contourf
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.contourf(Vorticity.T, 40, cmap='jet')
    cbar = fig.colorbar(ax.contourf(Vorticity.T, 40, cmap='jet'))
    cbar.ax.set_ylabel('Vorticity')
    ax.set_title("Vorticity")
    plt.savefig("Vorticity.png",dpi=300)


    ## Compute the Stream function on the interior points using the vorticity using SOR method for solving the Laplace equation laplacian(ψ) = -ω
    ## use omega = 1.5 and tolerance = 1e-6 for the SOR method
    Stream = np.zeros((N,N))

    tolerance = 1e-6
    max_iter = 10000
    iter = 0

    omega = 1.5
    residual = 100 ## entry criteria for the while loop
    while iter < max_iter and residual > tolerance:
        residual = 0
        StreamOld = Stream.copy()
        for i in range(1,N-1):
            for j in range(1,N-1):
                currValue = 0.25*(Stream[i+1,j] + Stream[i-1,j] + Stream[i,j+1] + Stream[i,j-1] + Vorticity[i,j])
                Stream[i,j] = (1-omega)*Stream[i,j] + omega*currValue
        residual = np.linalg.norm(Stream - StreamOld)
        
        ## print the residual every 100 iterations fixed width 10 and precision 6 using fstrings
        if iter % 100 == 0:
            print(f"iter: {iter:10d} residual: {residual:10.6f}")

        iter += 1
        
    ## plot the Stream function in a contour plot with contour levels 90 and colormap jet
    ## plot it using ax.contourf
    fig, ax = plt.subplots(1,1, figsize=(5,5))
    ax.contourf(Stream.T, 40, cmap='jet')
    ax.set_title("Stream Function")
    # Add colorbar
    cbar = fig.colorbar(ax.contourf(Stream.T, 40, cmap='jet'))
    cbar.ax.set_ylabel('Stream Function')
    plt.savefig("Stream.png",dpi=300)

# %%
