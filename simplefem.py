import numpy as np
from math import *
import drawMesh
import pygmsh
import gmsh

gmsh.initialize()

rect_width, rect_length  = 3.0, 10.0
resolution = 0.1

geom = pygmsh.geo.Geometry()


circle = geom.add_circle(
    [5,1.5,0],
    radius=0.5,
    mesh_size=resolution*0.5,
    make_surface=False
)


rect = geom.add_polygon(
    [
        [0.0, 0.0                ,0],
        [0.0, rect_width         ,0],
        [rect_length, rect_width ,0],
        [rect_length, 0.0        ,0],
    ],
    mesh_size=resolution , holes = [circle]
    
)

mesh = geom.generate_mesh(dim=2)

geom.__exit__()


# ## Node

class Node:
    
    def __init__(self, id, x, y):
        self.id = id
        self.x, self.y = x, y
        self.fx, self.fy = 0.0, 0.0
        self.rx, self.ry = 0.0 ,0.0
        self.dx, self.dy = None, None
        
    
    @property
    def dfix(self):
        if self.dx == 0.0 and self.dy == 0.0:
            return True
        else:
            return False
        
    @property
    def externalForce(self):
        if self.fx != 0.0 or self.fy != 0.0:
            return True
        else:
            return False
        
    def __eq__(self, obj):
        if (self.x == obj.x) and (self.y == obj.y):
            return True
        else:
            return False       
        
    def __str__(self):
        return str(self.__dict__)

# %% [markdown]
# ## Element

# %%
class Element:
    
    maxColorVal = -9.9e19
    minColorVal = 9.9e19
    colorFunc = lambda x: x
    
    def __init__(self, id, nodes):
        self.id = id
        self.nodes = self.orderCounterClock(nodes)
        self.stress = None
        self.strain = None
        self.colorVal = 0
        self.getArea()
    
    @property
    def getmaxColorVal(self):
        return Element.maxColorVal
    
    @property
    def getminColorVal(self):
        return Element.minColorVal
    
    @property
    def getcolorFunc(self):
        return Element.colorFunc

    def getde(self):
        de_ = []
        for n in self.nodes:
            de_.append(n.dx)
            de_.append(n.dy)
        self.de = np.array(de_)
        return self.de
    
    def getColor(self):
        
        try: x_ = float(self.colorVal - Element.minColorVal)/(Element.maxColorVal - Element.minColorVal)
        except ZeroDivisionError: x_ = 0.5 # cmax == cmin
        
        x = Element.colorFunc(x_)
        
        blue  = int(255* min((max((4*(0.75-x), 0.)), 1.)))
        red   = int(255* min((max((4*(x-0.25), 0.)), 1.)))
        green = int(255* min((max((4*fabs(x-0.5)-1., 0.)), 1.)))
        return (red, green, blue)
    
    def getArea(self):
        x1,y1 = self.nodes[0].x, self.nodes[0].y
        x2,y2 = self.nodes[1].x, self.nodes[1].y
        x3,y3 = self.nodes[2].x, self.nodes[2].y
        result = 0.5*((x2*y3 - x3*y2)-(x1*y3- x3*y1)+(x1*y2-x2*y1))
        if result == 0:
            result = 1e-20
        self.area = result
        return result
    
    def getBe(self):
        x1,y1 = self.nodes[0].x, self.nodes[0].y
        x2,y2 = self.nodes[1].x, self.nodes[1].y
        x3,y3 = self.nodes[2].x, self.nodes[2].y
        B = (0.5/self.area) * np.array([
            [(y2-y3) ,  0    , (y3-y1),  0   ,   (y1-y2),   0   ],
            [   0    , (x3-x2),  0    , (x1-x3),     0   ,(x2-x1)],
            [(x3-x2) , (y2-y3), (x1-x3), (y3-y1), (x2-x1) ,(y1-y2)],
        ], dtype=np.float64)
        self.Be = B
        return B
        
    def getKe(self, D):
        Bie = self.getBe()
        Ke = self.area* np.matmul(Bie.T, np.matmul(D, Bie))
        self.Ke = Ke
        return Ke
    
    def orderCounterClock(self, nodes):
        p1,p2,p3 = nodes[0], nodes[1], nodes[2]
        val = (p2.y - p1.y) * (p3.x - p2.x) - (p2.x - p1.x) * (p3.y - p2.y)
        nodes_ = nodes.copy()
        if val > 0:
            nodes[1] = nodes_[0]
            nodes[0] = nodes_[1]   
        
        assembly = []
        for n in nodes:
            assembly.append(int(n.id*2))
            assembly.append(int(n.id*2) +1)
        self.assembly = assembly
        
        return nodes
    
    def __str__(self):
        return str(self.id) + ': [ ' + ', '.join([str(node.id) for node in self.nodes]) + ' ]'
        

# %% [markdown]
# # Preprocessing

# %% [markdown]
# ## Extract mesh data

# %%
maxNode = 0
for cell in mesh.cells[1].data:
    for node in cell:
        if node > maxNode:
            maxNode = node

meshCells = mesh.cells[1].data - np.full(np.shape(mesh.cells[1].data), 1, dtype=np.uint64)
meshPoints = mesh.points[1:]

# %%
nodes = [Node(i, point[0], point[1]) for i, point in enumerate(meshPoints)]
elements = []

for i,cell in enumerate(meshCells):
    elements.append(
        Element(id=i, nodes=[nodes[i] for i in cell])
    )

# %% [markdown]
# ## Material properties

# %%
v = 0.28
E = 200.0e9

D = (E/(1-v**2)) * np.array([
    [1, v, 0],
    [v, 1, 0],
    [0, 0, (1-v)/2],
])

# %% [markdown]
# ## Boundary conditions and forces

# %%
for i, node in enumerate(nodes):
    if node.x == rect_length:           # At right side of the rectangle (x=10)
        node.fx = 1.0e3                 # Apply a downwards force of 1kN
    elif node.x == 0.0:                 # At left side of the rectangle (x=0)
        node.dx, node.dy = 0.0, 0.0     # Fix the displacement in x and y
        node.rx, node.ry = None, None   # Set the reaction forces as unknowns


# # Stiffness matrix assembly

# %%
def assemblyK(K, Ke, nodeAssembly):
    for i,t in enumerate(nodeAssembly):
        for j,s in enumerate(nodeAssembly):
            K[t][s] += Ke[i][j]

# %%
Nnodes = len(nodes)
K = np.zeros((Nnodes*2,Nnodes*2))
B_list = []

for e in elements:
    Ke = e.getKe(D)
    nodeAssembly = e.assembly
    assemblyK(K, Ke, nodeAssembly)

# %% [markdown]
# $\begin{equation}
# Kd = f + r
# \end{equation}$

# %%
f = np.zeros((int(2*Nnodes), 1))
d = np.full((int(2*Nnodes), 1), None)
r = np.full((int(2*Nnodes), 1), None)

rowsrk, rowsdk = [], []

for i,node in enumerate(nodes):
    ix,iy = int(i*2), int(i*2)+1
    
    f[ix], f[iy] = node.fx, node.fy
    d[ix], d[iy] = node.dx, node.dy
    r[ix], r[iy] = node.rx, node.ry
    
    if node.dx == None:
        rowsrk.append(ix)
    else:
        rowsdk.append(ix)
        
    if node.dy == None:
        rowsrk.append(iy)
    else:
        rowsdk.append(iy)

# %%
rowsrk = [i for i in range(len(d)) if d[i] == None]
rowsdk = [i for i in range(len(r)) if r[i] == None]

# %% [markdown]
# # Solver

# %%
KB = np.zeros((len(rowsrk),len(rowsrk)))
KA = np.zeros((len(rowsdk),len(rowsrk)))

fk = np.array([r[i] for i in rowsrk]) + np.array([f[i] for i in rowsrk])
dk = np.array([d[i] for i in rowsdk]) 

for i in range(np.shape(KB)[0]):
    for j in range(np.shape(KB)[1]):
        KB[i][j] = K [rowsrk[i]][rowsrk[j]]

for i in range(np.shape(KA)[0]):
    for j in range(np.shape(KA)[1]):
        KA[i][j] = K [rowsdk[i]][rowsrk[j]]


# %%
du = np.matmul(np.linalg.inv(KB), fk)
fu = np.matmul(KA,du)

# %% [markdown]
# # Postprocessing

# %%
d_total = d.copy()

for i, d_solve in zip(rowsrk, du):
    d_total[i] = d_solve

# %%
for i,n in enumerate(nodes):
    ix,iy = int(i*2), int(i*2)+1
    
    n.dx = d_total[ix][0]
    n.dy = d_total[iy][0]

# %% [markdown]
# ## Von-Mises Stress

# %%
def calculateVonMises(sx, sy, sxy):
    return sqrt(sx**2 + sy**2 + 3*(sxy**2) - sx*sy)

# %% [markdown]
# ## Mesh deformation and Coloring

# %%
def rgb(mag, cmin, cmax):
    
    try: x = float(mag-cmin)/(cmax-cmin)
    except ZeroDivisionError: x = 0.5 
    
    blue  = int(255* min((max((4*(0.75-x), 0.)), 1.)))
    red   = int(255* min((max((4*(x-0.25), 0.)), 1.)))
    green = int(255* min((max((4*fabs(x-0.5)-1., 0.)), 1.)))
    return (red, green, blue)

average = lambda x: (sum(x)/len(x))

# %%
maxd, mind = max(d_total)[0], min(d_total)[0]

Element.colorFunc = lambda x: x #exp(-x)

# %%
for i,element in enumerate(elements):
    
    de = element.getde()
    strain_e = np.matmul(element.Be,de)
    stress_e = np.matmul(D, strain_e)

    dx_avg = average([de[0], de[2], de[4]])
    dy_avg = average([de[1], de[3], de[5]])
    
    element.strain = strain_e
    element.stress = stress_e
    
    element.colorVal = calculateVonMises(element.stress[0], element.stress[1], element.stress[2])
    
    if element.colorVal > Element.maxColorVal:
        Element.maxColorVal = element.colorVal
    if element.colorVal < Element.minColorVal:
        Element.minColorVal = element.colorVal
        

# %%
render = drawMesh.MeshRender()
render.legend = True
render.autoScale = True
render.deform_scale = 1.0e5
render.legendDiscretize = 10
render.legendTitle = 'von-mises (Pa)'
render.drawElements(elements)


#%%