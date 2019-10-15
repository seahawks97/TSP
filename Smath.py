"""
Author: Steven Tucker
Purpose: Self-coded mathematical equations to be used (mostly for
CSC 340)
"""
import math
PI = 3.14159265359
TWO_PI = 2 * PI
e = 2.71828
IMAGINARY = 1j


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
    Finds the distance between two points with coordinates (x1, y1)
    and (x2, y2)
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
    Quadratic formula: ax^2 + bx + c = 0
    Gives back 1 solution, 2 real solutions, or 2 complex solutions.
    :param a: leading coefficient of quadratic eqn
    :type a: float or int
    :param b: coefficient of second term in quadratic eqn
    :type b: float or int
    :param c: coefficient of last term in quadratic eqn
    :type c: float or int
    :return x: real roots of equation
    :rtype x: int or float
    """
    discriminant = (b**2) - (4*a*c)

    if discriminant < 0:
        discriminant *= -1
        upper_sqrt = sqrt(discriminant) * IMAGINARY
    elif discriminant == 0:
        x = -1*b/(2*a)
        return x
    else:
        upper_sqrt = sqrt(discriminant)

    x1 = ((-1*b) + upper_sqrt) / (2*a)
    x2 = ((-1*b) - upper_sqrt) / (2*a)
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
    if type(a) is not int or type(b) is not int:
        raise TypeError('Input must be float type.')

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
    if type(a) is not int or type(b) is not int:
        raise TypeError('Input must be float type.')

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
    return sqrt(variance(nums, avg, exhaustive=False))


def factorial(n):
    """
    Tells the value of n!
    :param n: number to be factorial
    :type n: int
    :return: n!
    :rtype: int
    """
    if type(n) is not int:
        raise TypeError('n must be type int.')

    if n > 1:
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


def Fibonacci(n):
    """
    Finds the nth term in the Fibonacci Sequence by adding the two
    previous terms.
    :param n: index of value
    :type n: int
    :return: the nth term i the Fibonacci Sequence
    :rtype: int
    """
    if type(n) is not int:
        raise TypeError('n must be type int')
    if n < 0:
        raise ValueError('n must be positive')

    if n == 0:
        return 0
    elif n == 1:
        return 1

    return Fibonacci(n-1) + Fibonacci(n-2)


def cmag(z):
    """
    Magnitude of a complex number
    :param z: complex number to find magnitude
    :type z: complex
    :return: magnitude
    :rtype: float
    """
    return sqrt(z.real**2 + z.imag**2)


def cphase(z):
    """
    Gives the phase of the complex number
    :param z: a complex number
    :type z: complex
    :return: angle in 2d-coords, or phase [-pi, pi]
    :rtype: float (radians)
    """
    return math.atan2(z.imag, z.real)


def deg2rad(deg):
    """
    Converts degrees to radians
    :param deg: measurement of degrees
    :type deg: float or int
    :return: radians
    :rtype: float or int
    """
    return deg * PI / 180


def rad2deg(rad):
    """
    Converts radians to degrees
    :param rad: measurement of radians
    :type rad: float or int
    :return: degrees
    :rtype: float or int
    """
    return rad * 180 / PI
