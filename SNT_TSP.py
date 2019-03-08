from Matrices import *
import itertools
import Smath
import time

# Number of nodes/cities for problem
tsp_size = 14

def get_tsp_data():
    """
    Creates a usable Matrix of TSP nodes.
    [i][0] is x-coordinate of node
    [i][1] is y-coordinate of node
    :return DM: Node coordinates
    :rtype DM: Matrix
    """
    DM = Matrix(m=14, n=2)
    DM.rows[0][0], DM.rows[0][1] = 0.430796749, 0.660341257
    DM.rows[1][0], DM.rows[1][1] = 0.869607109, 0.145710154
    DM.rows[2][0], DM.rows[2][1] = 0.272249997, 0.281035268
    DM.rows[3][0], DM.rows[3][1] = 0.310050105, 0.717362826
    DM.rows[4][0], DM.rows[4][1] = 0.847481151, 0.505130257
    DM.rows[5][0], DM.rows[5][1] = 0.054921944, 0.597324847
    DM.rows[6][0], DM.rows[6][1] = 0.565507064, 0.578311901
    DM.rows[7][0], DM.rows[7][1] = 0.578311901, 0.636552793
    DM.rows[8][0], DM.rows[8][1] = 0.170565332, 0.160881561
    DM.rows[9][0], DM.rows[9][1] = 0.800726237, 0.384045138
    DM.rows[10][0], DM.rows[10][1] = 0.622627218, 0.957622127
    DM.rows[11][0], DM.rows[11][1] = 0.608021461, 0.736718151
    DM.rows[12][0], DM.rows[12][1] = 0.628737267, 0.622146623
    DM.rows[13][0], DM.rows[13][1] = 0.665929436, 0.720342005
    return DM


def create_distances(mat, size):
    """
    Creates a nested list of distance between nodes. The entry in the
    ith row and the jth column is the cost (distance) of traveling from
    city i to city j.
    :param mat: Node data coordinates
    :type mat: Matrix
    :return dist: Distance from city i to city j
    :rtype dist: Nested list
    """
    dist = [[Smath.distance(x1=mat.rows[i][0], y1=mat.rows[i][1],
                     x2=mat.rows[j][0], y2=mat.rows[j][1]) for j in range(size)] for i in range(size)]
    return dist

tsp_nodes = get_tsp_data()
distances = create_distances(tsp_nodes, tsp_size)
# print(distances)

def get_minmax(dist):
    """
    Gets the minimum and maximum distances from each city to a specific
    city. Determines a and b (min/max) in histogram.
    :param dist: distances between nodes
    :type dist: array (nested list)
    :return: min/max distances, and their corresponding cities
    :rtype: list
    """
    min, max = dist[0][1], dist[0][2]
    mini, minj, maxi, maxj = 0, 1, 0, 2
    for i in range(len(dist)):
        for j in range(len(dist[i])):
            cur = dist[i][j]
            if cur != 0:
                if cur > max:
                    max = cur
                    maxi, maxj = i, j
                if cur < min:
                    min = cur
                    mini, minj = i, j

    return [min, max, maxi, maxj, mini, minj]

minmax_vals = get_minmax(distances)
min_cost11 = minmax_vals[0]
max_cost11 = minmax_vals[1]
max_route11 = (minmax_vals[2], minmax_vals[3])# from 1st city to 2nd city
min_route11 = (minmax_vals[4], minmax_vals[5])# from 1st city to 2nd city
print('Min/Max vals (Jordan, ignore):', minmax_vals)


################################################################# Admin
# Original problem states that 14 cities are nodes and stops, but I
# used less to test how my algorithm behaves.
#
number_of_nodes = 10
number_of_stops = 6


############################################################# Histogram
################################### Still need to export these to Excel
# Absolute min/max round trip costs.
min_cost_rt = min_cost11 * number_of_nodes
max_cost_rt = max_cost11 * number_of_nodes
# Bin size depend on number of bins (user input) and abs min/max
n_bins = 128
bins = [0 for i in range(n_bins)]
cost_range = max_cost_rt - min_cost_rt
dx = cost_range / n_bins

########################################################## Permutations
# list of city labels
node_list = [int(i) for i in range(number_of_nodes)]

n = 0 # number of routes
sumx = 0 # sum of cost
sumxx = 0 # sum of cost squared
LC = 0 # Float, least cost
LCT = () # Tuple, route w LC
MC = 0 # Float, most cost
MCT = () # Tuple, route w MC

# This is it: The big loop
# Start timer
start = time.time()
for route in itertools.permutations(node_list, number_of_stops):
    # print('Route:', route)
    this_dist = 0

    # Find the cost of a round trip
    for stop in range(len(route)):
        if stop != len(route) - 1:      # if its not the last node
            this_dist += distances[route[stop]][route[stop+1]]
        else:       # Add cost from last item to first item
            this_dist += distances[route[stop]][route[0]]

    # print('Round trip distance', this_dist)

    # Initialize LeastCost and MostCost
    # if total distance is 0, ie before you actually visited any nodes
    # sets LC, MC to first cost
    # sets LCT, MCT to first route
    if sumx == 0:
        LC = this_dist
        LCT = route
        MC = this_dist
        MCT = route

    # Keep totals
    n += 1
    sumx += this_dist
    sumxx += this_dist ** 2

    # Update min/max
    if this_dist < LC:
        LC = this_dist
        LCT = route
    elif this_dist > MC:
        MC = this_dist
        MCT = route

    # Sort cost into bin for histogram of frequencies
    bins[int((this_dist - min_cost_rt)/dx)] += 1

# Stop timer
end = time.time()

# Output
print('\nSum:', sumx)
print('Sum squared x:', sumxx)
print('Number of routes:', n)
print('Average cost:', sumx/n)
print('Least cost:', LC)
print('Least cost trip:', LCT)
print('Most cost:', MC)
print('Most cost trip:', MCT)
print('Bins:', bins)
print('Time:', end-start, 'seconds')


