
pi = 3.14159265359
e = 2.71828

def sqrt(n):
    """
    Finds the square root of the number given.
    :param n: The number being square rooted
    :type n: float or int
    :returns: square-rooted number
    :rtype: float
    """
    return n ** 0.5


def distance(x1, y1, x2, y2):
    """
    Finds the distance between two points with coordinates (x1, y1) and (x2, y2)
    :param x1: x-coordinate of first point
    :type x1: float or int
    :param y1: y-coordinate of first point
    :type y1: float or int
    :param x2: x-coordinate of second point
    :type x2: float or int
    :param y2: y-coordinate of second point
    :type y2: float or int
    :return: distance between two points
    :rtype: float
    """
    return sqrt(((x2-x1)**2) + ((y2-y1)**2))


def quadEqnSolver(a, b, c):
    """
    Applies the quadratic formula to given coefficients.
    :param a: leading coefficient of quadratic eqn
    :type a: float or int
    :param b: coefficient of second term in quadratic eqn
    :type b: float or int
    :param c: coefficient of last term in quadratic eqn
    :type c: float or int
    :return x: real roots of equation
    :rtype x: int or float
    """
    # print('Follow the form of the equation ax^2 + bx + c = 0')
    # a = int(input('Enter the coefficient for a: '))
    # b = int(input('Enter the coefficient for b: '))
    # c = int(input('Enter the coefficient for c: '))
    d = (b**2) - (4*a*c)

    if d < 0:
        # print('Zero real solutions')
        return None
    elif d == 0:
        x = -1*b/(2*a)
        # print('One solution:', x)
        return x
    else:
        x1 = (-1*b + sqrt(d)) / (2*a)
        x2 = (-1*b - sqrt(d)) / (2*a)
        # print('Two solutions:', x1, ',', x2)
        return x1, x2


def LCM(a, b):
    """
    Finds the Least Common Multiple of two given numbers
    :param a: arbitrary first number
    :type a: int
    :param b: arbitrary second number
    :type b: int
    :return: least common multiple
    :rtype: int
    """
    return a / GCF(a, b) * b


def GCF(a, b):
    """
    Finds Greatest Common Factor of two given numbers
    :param a: arbitrary first number
    :type a: int
    :param b: arbitrary second number
    :type b: int
    :return: greatest common factor
    :rtype: int
    """
    if b > a:
        return GCF(b, a)

    if a % b == 0:
        return b

    return GCF(b, a % b)


def mean(nums):
    """
    Gets mean value of a list of numbers
    :param nums: contains numbers to be averaged
    :type nums: list
    :return: average of nums
    :rtype: float or int
    """
    counter = 0
    for i in nums:
        counter += i
    return counter / len(nums)


def variance(nums, avg, exhaustive=False):
    """
    Gives variance of input data
    :param nums: contains numbers to use
    :type nums: list
    :param avg: average of nums
    :type avg: float or int
    :param exhaustive: optimizes how data is retrieved in TSP
    :type exhaustive: bool
    :return: variance
    :rtype: float or int
    """
    if not exhaustive:
        total = 0
        for i in nums:
            total += (i - avg) ** 2
        return total / len(nums)

    else:
        t1, t2 = 0, 0
        for i in nums:
            t1 += i ** 2
            t2 += i
        return (t1 - ((t2 ** 2) / len(nums))) / len(nums)


def standard_deviation(nums, avg):
    """
    Gives standard deviation of input data
    :param nums: contains numbers to use
    :type nums: list
    :param avg: average of nums
    :return: size of one standard deviation
    :rtype; float or int
    """
    return sqrt(variance(nums, avg, exhaustive = False))

def factorial(n):
    """
    Tells the value of n!
    :param n: number to be factorial
    :type n: int
    :return: n!
    :rtype: int
    """
    if n >= 0 and type(n) is int:
        if n == 0:
            return 1
        else:
            return n * factorial(n-1)
    else:
        return 1


def nPr(n, r):
    """
    Permuter gives number of permutations of r numbers from a selection
    of n points. Order does matter.
    :param n: total number of data points available to use
    :type n: int
    :param r: number of points actually used
    :type r: int
    :return: number of permutations using r numbers of n points
    :rtype: int
    """
    return int(factorial(n)/factorial(n-r))


def nCr(n, r):
    """
    Chooser gives number of combinations of r numbers from a selection
    of n points. Order does not matter.
    :param n: total number of ddata points available to use
    :type n: int
    :param r: number of points actually used
    :type r: int
    :return: number of combinations using r numbers of n points
    :rtype: int
    """
    return int(factorial(n)/(factorial(r) * factorial(n-r)))

