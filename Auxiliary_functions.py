import numpy as np
from typing import Dict, List, Tuple, Optional
import scipy.sparse as sp

def build_grid(Nx: int, Ny: int, Lx: float, Ly: float) -> Tuple[np.ndarray, np.ndarray, float, float]:
    x = np.linspace(0.0, Lx, Nx)
    y = np.linspace(0.0, Ly, Ny)
    dx = x[1] - x[0]
    dy = y[1] - y[0]
    return x, y, dx, dy

def idx(i: int, j: int, Ny: int) -> int:
    '''
    Indexing (i,j) -> k, with i in x and j in y, shape (Nx,Ny).
    '''
    return i * Ny + j

def build_laplacian_periodic_x_neumann_y(Nx: int, Ny: int, dx: float, dy: float) -> sp.csr_matrix:
    '''
    Laplacian on a staggered grid:
    - Periodicity in x
    - Neumann condition in y.
    '''
    N = Nx * Ny
    rows, cols, data = [], [], []

    invdx2 = 1.0 / dx**2
    invdy2 = 1.0 / dy**2

    for i in range(Nx):
        iL = (i - 1) % Nx
        iR = (i + 1) % Nx
        for j in range(Ny):
            k = idx(i, j, Ny)

            # Coef central:
            cP = -2.0 * (invdx2 + invdy2)

            # x neighbors (periodic)
            kL = idx(iL, j, Ny)
            kR = idx(iR, j, Ny)

            # y neighbors (Neumann)
            if j == 0:
                # down: phi_{j-1} = phi_{j+1}
                kS = idx(i, j + 1, Ny)
                kN = idx(i, j + 1, Ny)
            elif j == Ny - 1:
                # Top: phi_{j+1} = phi_{j-1}
                kS = idx(i, j - 1, Ny)
                kN = idx(i, j - 1, Ny)
            else:
                kS = idx(i, j - 1, Ny)
                kN = idx(i, j + 1, Ny)

            # ensamble
            rows += [k, k, k, k, k]
            cols += [k, kL, kR, kS, kN]
            data += [cP, invdx2, invdx2, invdy2, invdy2]

    L = sp.csr_matrix((data, (rows, cols)), shape=(N, N))
    return L

def divergence(u: np.ndarray, v: np.ndarray, dx: float, dy: float) -> np.ndarray:
    '''
    Centered nabla*u (periodicity in x is implicit if you wrap in i).
    '''
    Nx, Ny = u.shape
    div = np.zeros_like(u)

    for i in range(Nx):
        iL = (i - 1) % Nx
        iR = (i + 1) % Nx
        for j in range(1, Ny - 1):
            dudx = (u[iR, j] - u[iL, j]) / (2.0 * dx)
            dvdy = (v[i, j + 1] - v[i, j - 1]) / (2.0 * dy)
            div[i, j] = dudx + dvdy

    return div

def vorticity(u: np.ndarray, v: np.ndarray, dx: float, dy: float) -> np.ndarray:
    '''
    w = dv/dx - du/dy (centered).
    '''
    Nx, Ny = u.shape
    w = np.zeros_like(u)
    for i in range(Nx):
        iL = (i - 1) % Nx
        iR = (i + 1) % Nx
        for j in range(1, Ny - 1):
            dvdx = (v[iR, j] - v[iL, j]) / (2.0 * dx)
            dudy = (u[i, j + 1] - u[i, j - 1]) / (2.0 * dy)
            w[i, j] = dvdx - dudy
    return w

