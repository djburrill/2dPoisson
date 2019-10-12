# poissonRelax.py

# Imports
import numpy as np
import matplotlib.pyplot as plt

# Functions
def gaussian(x,y,cx,cy):
    '''
    Return Gaussian value at (x,y) while centered at (cx,cy)
    '''

    # Variables
    A = 1.0
    sigX = 2.0
    sigY = 2.0

    # Calculate terms
    term1 = ((x-cx)*(x-cx))/(2.0*sigX*sigX)
    term2 = ((y-cy)*(y-cy))/(2.0*sigY*sigY)

    # Return Gaussian
    return A*np.exp(-(term1+term2))

def poissonRelax():
    '''
    Perform relaxation method on Poisson equation
    '''

    # Variables
    sqExtent = 10
    numGridSide = 500
    gridSize = 1.0*sqExtent/numGridSide
    tol = 0.001                         # Define change tolerance
    diff = 2*tol                        # Difference parameter

    # Function parameters
    centers = [[3.0,8.0],[7.0,2.0]]
    funcVal = 0.0

    # Initialize grid
    grid = np.zeros((numGridSide,numGridSide))

    # Relaxation loop
    while(diff > tol):
        #print 'Diff: ' + str(diff)
        diff = 0.0
        # Loop over grid points
        for index1 in range(numGridSide):
            for index2 in range(numGridSide):
                # Store initial value
                origVal = grid[index2,index1]

                # Calculate X and Y
                X = index1*gridSize
                Y = index2*gridSize

                # Calculate function at position
                funcVal = 0.0

                for centerPos in centers:
                    funcVal += gaussian(X,Y,centerPos[0],centerPos[1])

                # Calculate new grid value
                # Bottom left corner
                if (index2 == 0 and index1 == 0):
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2,index1+1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Bottom right corner
                elif (index2 == 0 and index1 == numGridSide-1):
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Top right corner
                elif (index2 == numGridSide-1 and index1 == numGridSide-1):
                    grid[index2,index1] = (grid[index2-1,index1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Top left corner
                elif (index2 == numGridSide-1 and index1 == 0):
                    grid[index2,index1] = (grid[index2-1,index1] +
                                           grid[index2,index1+1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Left side
                elif (index1 == 0):
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2-1,index1] +
                                           grid[index2,index1+1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Right side
                elif (index1 == numGridSide-1):
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2-1,index1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Top side
                elif (index2 == numGridSide-1):
                    grid[index2,index1] = (grid[index2-1,index1] +
                                           grid[index2,index1+1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Bottom side
                elif (index2 == 0):
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2,index1+1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25
                # Center
                else:
                    grid[index2,index1] = (grid[index2+1,index1] +
                                           grid[index2-1,index1] +
                                           grid[index2,index1+1] +
                                           grid[index2,index1-1] -
                                           gridSize*gridSize*funcVal)*0.25

                # Check difference
                if (np.fabs(origVal-grid[index2,index1]) > diff):
                    diff = np.fabs(origVal-grid[index2,index1])

    # Plot results
    plt.imshow(grid, cmap='PuBu')
    plt.title('Poisson Equation')
    plt.xticks([])
    plt.yticks([])
    plt.show()

# Main
if (__name__ == '__main__'):
    poissonRelax()
