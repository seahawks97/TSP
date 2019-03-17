import itertools
import Smath
import time
import random

# Number of nodes/cities for problem
tsp_size = 14

def main():
    # Setup
    tsp_nodes = get_tsp_data()
    distances = create_distances(tsp_nodes, tsp_size)
    # print(distances)

    minmax_vals = get_minmax(distances)
    min_cost11 = minmax_vals[0]
    max_cost11 = minmax_vals[1]
    # max_route11 = (minmax_vals[2], minmax_vals[3])  # from 1st city to 2nd city
    # min_route11 = (minmax_vals[4], minmax_vals[5])  # from 1st city to 2nd city
    # print('Min/Max vals:', minmax_vals)

    # Admin
    number_of_nodes = int(input('How many cities are there (1 <= x <= 14): '))
    number_of_stops = int(input('How many to visit (2 <= x <= above): '))

    # Histogram
    # Absolute min/max round trip costs.
    min_cost_rt = min_cost11 * number_of_nodes
    max_cost_rt = max_cost11 * number_of_nodes
    # Bin size depend on number of bins (user input) and abs min/max
    n_bins = int(input('How many bins for histogram ( x >= 100): '))
    bins = [0 for i in range(n_bins)]
    cost_range = max_cost_rt - min_cost_rt
    dx = cost_range / n_bins

    # Set stats up
    # list of city labels
    node_list = [int(i) for i in range(number_of_nodes)]
    sumx = 0  # sum of cost
    sumxx = 0  # sum of cost squared
    # LC = 0  # Float, least cost
    # LCT = ()  # Tuple, route w LC
    # MC = 0  # Float, most cost
    # MCT = ()  # Tuple, route w MC

    # Big loop
    selection = input('1. Exhaustive Search\n'
                      '2. Random Search\n'
                      '3. Genetic Algorithm\n'
                      '4. Simulated Annealing\n'
                      'Select an option above by typing the number and pressing Enter: ')
    if selection == '1':
        exhaustive_search(node_list, number_of_stops, distances, sumx, sumxx, bins, min_cost_rt, dx, number_of_nodes)
    elif selection == '2':
        randomized(distances, sumx, sumxx, bins, min_cost_rt, dx)
    elif selection == '3':
        pass
    elif selection == '4':
        pass

    # Stop exiting
    input('Press Enter to quit: ')


def get_tsp_data():
    """
    Creates a usable array of TSP nodes.
    [i][0] is x-coordinate of node
    [i][1] is y-coordinate of node
    :return DM: Node coordinates
    :rtype DM: Array
    """
    DM = [[0.430796749, 0.660341257],
          [0.869607109, 0.145710154],
          [0.272249997, 0.281035268],
          [0.310050105, 0.717362826],
          [0.847481151, 0.505130257],
          [0.054921944, 0.597324847],
          [0.565507064, 0.578311901],
          [0.578311901, 0.636552793],
          [0.170565332, 0.160881561],
          [0.800726237, 0.384045138],
          [0.622627218, 0.957622127],
          [0.608021461, 0.736718151],
          [0.628737267, 0.622146623],
          [0.665929436, 0.720342005]]
    return DM


def create_distances(mat, size):
    """
    Creates a nested list of distance between nodes. The entry in the
    ith row and the jth column is the cost (distance) of traveling from
    city i to city j.
    :param mat: Node data coordinates
    :type mat: Array
    :param size: Number of nodes
    :type size: Int
    :return dist: Distance from city i to city j
    :rtype dist: Array
    """
    dist = [[Smath.distance(x1=mat[i][0], y1=mat[i][1],
                     x2=mat[j][0], y2=mat[j][1]) for j in range(size)] for i in range(size)]
    return dist


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


def exhaustive_search(node_list, number_of_stops, distances, sumx, sumxx, bins, min_cost_rt, dx, number_of_nodes):
    """
    Exhaustive search for TSP. Simulates all possible combinations of cities chosen.
    :param node_list:
    :param number_of_stops:
    :param distances:
    :param sumx:
    :param sumxx:
    :param bins:
    :param min_cost_rt:
    :param dx:
    :param number_of_nodes:
    :return:
    """

    # Initialize LC, MC
    LC = 100
    MC = 0

    # This is it: The big loop
    # Start timer
    start = time.time()

    for route in itertools.permutations(node_list, number_of_stops):

        if route[0] < route[-1]:
            # print('Route:', route)
            this_dist = 0

            # Find the cost of a round trip
            for stop in range(len(route)):
                if stop != len(route) - 1:      # if its not the last node
                    this_dist += distances[route[stop]][route[stop+1]]
                else:       # Add cost from last item to first item
                    this_dist += distances[route[stop]][route[0]]

            # print('Round trip distance', this_dist)

            # Keep totals
            sumx += this_dist
            sumxx += this_dist ** 2

            # Update min/max
            if this_dist < LC:
                LC = this_dist
                LCT = route
                print('New low:', LCT, 'Distance:', LC)
            elif this_dist > MC:
                MC = this_dist
                MCT = route
                print('New high:', MCT, 'Distance:', MC)

            # Sort cost into bin for histogram of frequencies
            bins[int((this_dist - min_cost_rt)/dx)] += 1

            # print('Route:', route)

    # Stop timer
    end = time.time()

    print('\nSum:', sumx)
    print('Sum squared x:', sumxx)
    n = Smath.nPr(number_of_nodes, number_of_stops)
    print('Number of routes:', n)
    print('Average cost:', sumx / n)
    stdev = Smath.sqrt((sumxx - ((sumx ** 2) / n)) / n)
    print('Standard deviation:', stdev)
    print('Least cost:', LC)
    print('Least cost trip:', LCT)
    print('Most cost:', MC)
    print('Most cost trip:', MCT)
    print('Bins:', bins)
    print('Time:', end - start, 'seconds')


def randomized(distances, sumx, sumxx, bins, min_cost_rt, dx):
    """
    Generates bigN number of randomized trip orders. Checks the
    :param distances:
    :param sumx:
    :param sumxx:
    :param bins:
    :param min_cost_rt:
    :param dx:
    :return:
    """

    route = [i for i in range(14)]
    bigN = 1000
    LC = 100
    MC = 0

    # Start timer
    start = time.time()

    for j in range(bigN):
        this_dist = 0
        route = shuffle(route)
        # print('Route:', route)

        # Get the distance for the route
        for stop in range(len(route)):
            if stop != len(route) - 1:  # if its not the last node
                this_dist += distances[route[stop]][route[stop + 1]]
            else:  # Add cost from last item to first item
                this_dist += distances[route[stop]][route[0]]
        # print('Distance:', this_dist)

        # Keep totals
        sumx += this_dist
        sumxx += this_dist ** 2

        # Update min/max
        if this_dist < LC:
            LC = this_dist
            LCT = route
            print('New low:', LCT, 'Distance:', LC)
        elif this_dist > MC:
            MC = this_dist
            MCT = route
            print('New high:', MCT, 'Distance:', MC)

        # Sort cost into bin for histogram of frequencies
        bins[int((this_dist - min_cost_rt) / dx)] += 1

    # Stop timer
    stop = time.time()
    print(route, this_dist)

    # Output
    print('\nSum:', sumx)
    print('Sum squared x:', sumxx)
    print('Number of routes:', bigN)
    print('Average cost:', sumx / bigN)
    stdev = Smath.sqrt((sumxx - ((sumx ** 2) / bigN)) / bigN)
    print('Standard deviation:', stdev)
    print('Least cost:', LC)
    print('Least cost trip:', LCT)
    print('Most cost:', MC)
    print('Most cost trip:', MCT)
    print('Bins:', bins)
    print('Time:', stop - start, 'seconds')


def shuffle(idxs):
    """
    Performs random swaps n times with random seeded generator
    :param idxs: Set of items to choose from
    :type idxs: list, tuple, or dict
    :return: modified idxs
    :rtype modified idxs: type(idxs)
    """
    n = 32
    for k in range(n):
        random.seed()
        a, b = random.randint(0, 13), random.randint(0, 13)
        swap(idxs, a, b)
        # print(idxs)
    return idxs


def swap(mylist, a, b):
    """
    Swaps two values
    :param mylist: order of things to choose from
    :type mylist: list, tuple, or dict
    :param a: first index
    :type a: list/tuple: int; dict: key
    :param b: second index
    :type b: list/tuple: int; dict: key
    :return: none; edits items in place
    """
    temp = mylist[a]
    mylist[a] = mylist[b]
    mylist[b] = temp



main()