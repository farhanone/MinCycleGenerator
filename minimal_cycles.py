import numpy as np
import pandas as pd
import scipy
import vtk
import gudhi
import sys

def adj_loc_radi_to_csv(adj_file_name, solved_file_name, output_file_name="results/appended_adjancy_file.csv"):
    """
    This function takes spatial locations of nodes in spatially embedded graph
    along with correponding adjancy matrix. It loads the adjancy matrix into
    pandas dataframe and appends this dataframe with the columns containing
    spatial information of nodes.
    
    param: adj_file_name: Adjacancy file (must be readable from numpy.loadtxt())
    param: solved_file_name: file with spatial information of nodes
    param  output_file_name: path to store ourput file in csv formate
    """
    
    def get_inner_value(x):
        return x[0][0]
    
    # Load Adjacancy Matrix
    adj_matrix = np.loadtxt(adj_file_name, delimiter=',')
    # print(np.max(adj_matrix))
    df = pd.DataFrame(adj_matrix, columns = ["{:04n}".format(i) for i in range(len(adj_matrix))])
    
    # We invert the field, so that when using Gudhi library for filteration, this it adds edges with highest weights first.
    # Since the Gudhi do filtration from 0 to infinity, so by inverting it do our filteration as Max to 0.
    mul_inv =  max(df.max()) - df 
        
    mul_inv.values[[np.arange(mul_inv.shape[0])], [np.arange(mul_inv.shape[0])]] = 0 # make major diagonal zeros
    
    # Load location information from soved file
    mat = scipy.io.loadmat(solved_file_name)
    df1 = pd.DataFrame(mat['pres'][0]) # to pandas frame
    
    # Append new information (location, radii and id ) to the adjacancy info dataframe
    mul_inv['x'] = df1.x.apply(get_inner_value)
    mul_inv['y'] = df1.y.apply(get_inner_value)
    mul_inv['z'] = 0
    mul_inv['id'] = df1.id.apply(get_inner_value)
    mul_inv['radius'] = df1.r.apply(get_inner_value)
    mul_inv['radius_m'] = df1.rm.apply(get_inner_value)
    
    mul_inv.to_csv(output_file_name, index=False)


def get_VR(data, out_file='results/VR.vtp'):
    # This function create a Vietoris Rips using Gudhi library given a adjancy matrix
    # data: is a pandas dataframe, The distance matrix with x,y,z and id value extended on column
    filtration_ = []
    simplex_ = []
    
    distance_matrix = data.iloc[:, :-6].to_numpy()
    lookup_id = data.iloc[:, -6:-3].to_numpy()
    rips = gudhi.RipsComplex(distance_matrix = distance_matrix)
    simplex_tree = rips.create_simplex_tree(max_dimension=2)
    simplex_tree.prune_above_filtration(np.max(distance_matrix) - 1.00e-12)
    
    points = vtk.vtkPoints()
    edges =  vtk.vtkCellArray()
    triangles = vtk.vtkCellArray()
    edges_and_triangles = vtk.vtkCellArray()
    thre_arr = vtk.vtkFloatArray()
    thre_arr.SetName('filteration')
    tri_thre = []
    edge_thre = []
    for i in range(len(lookup_id)):
        points.InsertNextPoint(lookup_id[i])
    count = 0
    
    for simplex_with_filtration in simplex_tree.get_filtration():
        simplex = simplex_with_filtration[0]
        if len(simplex) == 2:
            edge = vtk.vtkLine()
            edge.GetPointIds().SetId(0, simplex[0])
            edge.GetPointIds().SetId(1, simplex[1])
            # edges_and_triangles.InsertNextCell ( edge );
            edges.InsertNextCell ( edge )
            # thre.InsertNextTuple1(simplex[1])
            edge_thre.append(simplex_with_filtration[1])
        if len(simplex) == 3:
            triangle = vtk.vtkTriangle()
            triangle.GetPointIds().SetId(0, simplex[0])
            triangle.GetPointIds().SetId(1, simplex[1])
            triangle.GetPointIds().SetId(2, simplex[2])
            # edges_and_triangles.InsertNextCell(triangle)
            triangles.InsertNextCell(triangle)
            # thre.InsertNextTuple1(simplex[1])
            tri_thre.append(simplex_with_filtration[1])
        # print(count)
        count = count+1

    thre = edge_thre + tri_thre
    for i in thre:
        thre_arr.InsertNextTuple1(i)


    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(edges)
    polydata.SetPolys(triangles)
    polydata.GetCellData().AddArray(thre_arr)



    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(out_file)
    w.SetInputData(polydata)
    w.Write()
    return simplex_tree, lookup_id



def add_edge(adj, src, dest):
    # Python3 code for printing shortest path between two vertices of unweighted graph
    # utility function to form edge between two vertices source and dest
    adj[src].append(dest)
    adj[dest].append(src)


def BFS(adj, src, dest, v, pred, dist):
    # A modified version of BFS that stores predecessor of each vertex in array p
    # and its distance from source in array d

    # a queue to maintain queue of vertices whose
    # adjacency list is to be scanned as per normal
    # DFS algorithm
    queue = []

    # boolean array visited[] which stores the
    # information whether ith vertex is reached
    # at least once in the Breadth first search
    visited = [False for i in range(v)]

    # initially all vertices are unvisited
    # so v[i] for all i is false
    # and as no path is yet constructed
    # dist[i] for all i set to infinity
    for i in range(v):

        dist[i] = sys.float_info.max
        pred[i] = -1

    # now source is first to be visited and
    # distance from source to itself should be 0
    visited[src] = True
    dist[src] = 0
    queue.append(src)

    # standard BFS algorithm
    while (len(queue) != 0):
        u = queue[0]
        queue.pop(0)
        for i in range(len(adj[u])):

            if (visited[adj[u][i]] == False):
                visited[adj[u][i]] = True
                dist[adj[u][i]] = dist[u] + 1
                pred[adj[u][i]] = u
                queue.append(adj[u][i])

                # We stop BFS when we find
                # destination.
                if (adj[u][i] == dest):
                    return True
    return False


def printShortestDistance(adj, s, dest, v):
    # Utility function to print the shortest distance
    # between source node and destination node

    # predecessor[i] array stores predecessor of i
    # and distance array stores distance of i from s
    pred=[0 for i in range(v)]
    dist=[0 for i in range(v)]

    if (BFS(adj, s, dest, v, pred, dist) == False):
        print("Given source and destination are not connected")
        return 0
    
    # vector path stores the shortest path
    path = []
    crawl = dest
    path.append(crawl)
    while (pred[crawl] != -1):
        path.append(pred[crawl])
        crawl = pred[crawl]

    # distance from source is in distance array
    # print("Shortest path length is : " + str(dist[dest]), end = '')

    # printing path from source to destination
    # print("\nPath is : : ")
    r=[]
    for i in range(len(path)-1, -1, -1):
        r.append(path[i])
    #print(path[i], end=' ')
    return  r

def get_cycles(simplex_tree):
    # Get cycles birth and death value
    cycles = []
    # min_persistence is the persistence threshold.
    # 0 means we are gonna ignore cycles with 0 persistence value.
    # To show all cycles for all persistence values assign -1.
    birth_death = simplex_tree.persistence(min_persistence=0.0) 
    for i in range(len(birth_death)):
        if birth_death[i][0] == 1:
            cycles.append(birth_death[i][1])
    print("Total number of cycles in VR : ", len(cycles))
    cycles.sort()
    return cycles

def get_simplex_indices(filtration):
    # Given cycle birth filteration value, get the newly added simplices.
    # (get simplex which are added at filteration value at birth cycle[0][0])
    # simplex_ = [list(x)[0] for x in simplex_tree.get_filtration()]
    # filteration_ = [list(x)[-1] for x in simplex_tree.get_filtration()]

    unique_filteration = set(filtration)
    return {x: [i for i, value in enumerate(filtration) if value == x] for x in unique_filteration}
    # simplex_[ind[4.0][1]]

# This function compute VR for all the filteration values untill given threshold
# It will store all the simplicies untill given thresdhold
def get_simplices_at_th(simplex_tree, threshold = 8):
    simplices = []
    filtration = []
    for x in simplex_tree.get_filtration():
        if x[-1] <= threshold and len(x[0]) > 1: 
            simplices.append(list(x)[0])
            filtration.append(list(x)[-1])
    unique_filteration = set(filtration)
    # print(unique_filteration)
    ind = {x: [i for i, value in enumerate(filtration) if value == x] for x in unique_filteration}
    return simplices, ind


# Get Adjacancy list-of-list for Breath-First-Search(BFS) algorithms
# algorithm taken from https://www.geeksforgeeks.org/shortest-path-unweighted-graph/
def get_adj_list(simplices, num_pts):
    adj_list = [[] for i in range(num_pts)];
    for i in simplices:
        if len(i)==2:
            add_edge(adj_list, i[0], i[1])
    return adj_list

def extract_cycle(simplex_tree, num_pts):
    # H1
    cycles = get_cycles(simplex_tree)

    # Below code will first identify the edges that are added at certain filteration values.
    # identifying that edge, we will remove that edge from adjancy list and then do BFS being ends points of edge as source and destination for BFS
    # This will create a dictionary of edges with their key as circle ID and 'birth' as their cycle.
    cycles_dict = {}
    persistence_value = []
    num_triangles=0
    for cycid in range(len(cycles)): # cyc[0] is the birth and cyc[1] is death of an cycle, cyc is an loop go through all the cycles
        # print(cycid)
        
        simplices, ind = get_simplices_at_th(simplex_tree, threshold = cycles[cycid][0])
        adj = get_adj_list(simplices, num_pts)
        # ind = get_simplex_indices(filtration)
        
        valid_cycle_birth = []
        valid_cycle_death = []
        
        for i in ind[cycles[cycid][0]]: # This loop is go through all the indexes in simplex list that are born at given filter value(cyc[0]
            # # Check if edges add at birth are also part of the death simplex which makes tham valid cycle
            # if cycles[cycid][1] == np.inf:
            # print("from ", adj[simplex_[i][1]]," remove ", simplex_[i][0])
            # print("from ", adj[simplex_[i][0]]," remove ", simplex_[i][1])
            adj[simplices[i][1]].remove(simplices[i][0])
            adj[simplices[i][0]].remove(simplices[i][1])
            # print("src ",simplex_[i][0], "dest", simplex_[i][1])

            r = printShortestDistance(adj, simplices[i][0], simplices[i][1],num_pts)
            if r == 0:
                adj[simplices[i][1]].append(simplices[i][0])
                adj[simplices[i][0]].append(simplices[i][1])
                continue
            if len(r) == 3: # sometime we got triangles and BFS takes them as cycle. so we remove one edge of triangle and run again.This will gives other representative cycle.
                # print('Triangle_cycles : ', cycid)
                num_triangles = num_triangles + 1
                adj[r[1]].remove(r[2])
                adj[r[2]].remove(r[1])

                r1 = printShortestDistance(adj, simplices[i][0], simplices[i][1],num_pts)# its r1 because the index is changes

                adj[r[2]].append(r[1])
                adj[r[1]].append(r[2])

                r = r1

            adj[simplices[i][1]].append(simplices[i][0])
            adj[simplices[i][0]].append(simplices[i][1])


            valid_cycle_birth.append(r)
            
            if cycles[cycid][1] == np.inf:
                valid_cycle_death.append(np.inf)
            else:
                # print(cycles[cycid][1])
                # # break
                # for j in ind[cycles[cycid][1]]: 
                #     if all([elem in simplices[j] for elem in  simplices[i]]):
                valid_cycle_death.append(simplices[i])

        if r==0:
            continue
        cycles_dict[cycid] = {'birth' : valid_cycle_birth, 'death' : valid_cycle_death, "persistance_value" : [cycles[cycid][0], cycles[cycid][1]]}
    print("Number of Triangles : ", num_triangles)
    return cycles_dict

# Given dict of cysles store them into VTP file. dict have all the points that makes the cycles.
def write_cycles_vtp(cycles, lookup_id, output_name):
    points = vtk.vtkPoints()
    for i in range(len(lookup_id)):
        points.InsertNextPoint(lookup_id[i])

    edges =  vtk.vtkCellArray()

    cycleId = vtk.vtkIntArray()
    cycleId.SetName("cycle_id")

    f = vtk.vtkFloatArray()
    f.SetName("birth_filter_thr")

    d = vtk.vtkFloatArray()
    d.SetName("death_filter_thr")

    cicle_size = vtk.vtkIntArray()
    cicle_size.SetName("cycle_size")
    
    responsible_edge = vtk.vtkIntArray()
    responsible_edge.SetName("responsible_edge")

    for key, value in cycles.items():
        # print(key)
        for i in range(len(value['birth'][0])-1):

            edge = vtk.vtkLine()
            edge.GetPointIds().SetId(0, value['birth'][0][i])
            edge.GetPointIds().SetId(1, value['birth'][0][i+1])
            edges.InsertNextCell( edge );
            cycleId.InsertNextTuple1(key)
            cicle_size.InsertNextTuple1(len(value['birth'][0]))
            f.InsertNextTuple1(value['persistance_value'][0])
            d.InsertNextTuple1(value['persistance_value'][1])
            responsible_edge.InsertNextTuple1(-1)


        edge = vtk.vtkLine()
        edge.GetPointIds().SetId(0, value['birth'][0][0])
        edge.GetPointIds().SetId(1, value['birth'][0][-1])
        edges.InsertNextCell( edge )
        cycleId.InsertNextTuple1(key)
        f.InsertNextTuple1(value['persistance_value'][0])
        d.InsertNextTuple1(value['persistance_value'][1])
        cicle_size.InsertNextTuple1(len(value['birth'][0]))
        responsible_edge.InsertNextTuple1(key)



    polydata = vtk.vtkPolyData()
    polydata.SetPoints(points)
    polydata.SetLines(edges)
    polydata.GetCellData().AddArray(cycleId)
    polydata.GetCellData().AddArray(cicle_size)
    polydata.GetCellData().AddArray(f)
    polydata.GetCellData().AddArray(d)
    polydata.GetCellData().AddArray(responsible_edge)
    w = vtk.vtkXMLPolyDataWriter()
    w.SetFileName(output_name)
    w.SetInputData(polydata)
    w.Write()


def get_cc(simplex_tree):
    # H0 on edges
    cc = []
    birth_death = simplex_tree.persistence(min_persistence=-1) 
    # min_persistence is the persistence threshold.
    # 0 means we are gonna ignore cycles with 0 persistence value.
    # To show all cycles for all persistence values assign -1.
    for i in range(len(birth_death)):
        if birth_death[i][0] == 0:
            cc.append(birth_death[i][1])
    print("Total number of cycles in VR : ", len(cc))
    cc.sort()
    return cc

# Get Adjacancy list-of-list for Breath-First-Search(BFS) algorithms
# algorithm taken from https://www.geeksforgeeks.org/shortest-path-unweighted-graph/
def get_adj_list(simplices, num_pts):
    adj_list = [[] for i in range(num_pts)];
    for i in simplices:
        if len(i)==2:
            add_edge(adj_list, i[0], i[1])
    return adj_list

def extract_cc(simplex_tree, num_pts):
    
    cc = get_cc(simplex_tree)
    # Below code will first identify the edges that are added at certain filteration values.
    # identifying that edge, we will remove that edge from adjancy list and then do BFS being ends points of edge as source and destination for BFS
    # This will create a dictionary of edges with their key as circle ID and 'birth' as their cycle.
    cc_dict = {}
    persistence_value = []
    num_triangles=0
    for ccid in range(len(cc)): # cyc[0] is the birth and cyc[1] is death of an cycle, cyc is an loop go through all the cycles
        # print(cycid)
        
        simplices, ind = get_simplices_at_th(simplex_tree, threshold = cc[ccid][0])
        adj = get_adj_list(simplices, num_pts)
        # ind = get_simplex_indices(filtration)
        
        valid_cc_birth = []
        valid_cc_death = []
        
        for i in ind[cc[ccid][0]]: # This loop is go through all the indexes in simplex list that are born at given filter value(cyc[0]
            # # Check if edges add at birth are also part of the death simplex which makes tham valid cycle
            # if cycles[cycid][1] == np.inf:
            # print("from ", adj[simplex_[i][1]]," remove ", simplex_[i][0])
            # print("from ", adj[simplex_[i][0]]," remove ", simplex_[i][1])
            adj[simplices[i][1]].remove(simplices[i][0])
            adj[simplices[i][0]].remove(simplices[i][1])
            # print("src ",simplex_[i][0], "dest", simplex_[i][1])

            r = printShortestDistance(adj, simplices[i][0], simplices[i][1],num_pts)
            if r == 0:
                adj[simplices[i][1]].append(simplices[i][0])
                adj[simplices[i][0]].append(simplices[i][1])
                continue
            if len(r) == 3: # sometime we got triangles and BFS takes them as cycle.
                # so we remove one edge of triangle and run again.This will gives other representative cycle.
                # print('Triangle_cycles : ', cycid)
                num_triangles = num_triangles + 1
                adj[r[1]].remove(r[2])
                adj[r[2]].remove(r[1])

                r1 = printShortestDistance(adj, simplices[i][0], simplices[i][1],num_pts)# its r1 because the index is changes

                adj[r[2]].append(r[1])
                adj[r[1]].append(r[2])

                r = r1

            adj[simplices[i][1]].append(simplices[i][0])
            adj[simplices[i][0]].append(simplices[i][1])


            valid_cc_birth.append(r)
            
            if cc[ccid][1] == np.inf:
                valid_cc_death.append(np.inf)
            else:
                
                # print(cycles[cycid][1])
                # # break
                # for j in ind[cycles[cycid][1]]: 
                #     if all([elem in simplices[j] for elem in  simplices[i]]):
                valid_cc_death.append(simplices[i])
                        
        if r==0:
            continue
            
        cc_dict[ccid] = {'birth' : valid_cc_birth, 'death' : valid_cc_death, "persistance_value" : [cc[ccid][0], cc[ccid][1]]}
    print("Number of Triangles : ", num_triangles)
    return cc_dict