# Steven Tucker
# CSC 340-001
# Toolbox: Matrix Kit

import Smath

class MatrixError(Exception):
    """Error raised for errors within Matrix operations."""
    pass

class Matrix():
    def __init__(self, m, n, init=True, identity=False):
        """
        Matrix class constructor
        :param m: number of rows
        :type m: int
        :param n: number of columns
        :type n: int
        :param init: initialize 0s?
        :type init: bool
        :param identity: initialize 1s on diagonal
        :type identity: bool
        """
        if init:
            self.rows = [[0]*n for x in range(m)]
            # self.columns = [[0]*m for y in range(n)]
        else:
            self.rows = []
            # self.columns = []

        self.m = m
        self.n = n

        if identity and self.isSquareMat():
            index = 0
            for row in self.rows:
                row[index] = 1
                index += 1

    def __str__(self):
        """
        Displays values of the matrix in their position
        :return s: numerical layout
        :rtype s: string
        """
        s = '\n'.join([' '.join([str(item) for item in row]) for row in self.rows])
        return s

    def __repr__(self):
        """
        Prints a more descriptive output of the matrix
        :return rep: numeric representation of matrix
        :rtype rep: string
        """
        s = str(self.rows)
        rep = "%s \nRows: %s, Columns: %s\n" % (s, self.mrow(), self.ncol())
        return rep

    def __setitem__(self, row, value):
        """
        Set a value to a given position in place
        :param row: position to edit
        :type row: int
        :param value: value to be set
        :type value: float or int
        :return: none
        """
        if type(row) is int:
            if type(value) is int or type(value) is float or type(value) is list:
                self.rows[row] = value
            else:
                raise MatrixError('Value must be type int or float')
        else:
            raise MatrixError('Row must be type int')

    def __getitem__(self, item):
        """
        Gets the value of a position
        :param item: item desired
        :type item: float or int
        :return: value in matrix
        :rtype: float or int
        """
        if type(item) is int:
            return self.rows[item]
        else:
            raise MatrixError('Unable to retrieve Matrix[item]; item must be int type')

    def mrow(self):
        """
        Number of rows
        :return: number of rows
        :rtype: int
        """
        return self.m

    def ncol(self):
        """
        Number of columns
        :return: number of columns
        :rtype: int
        """
        return self.n

    def getSize(self):
        """
        Dimensions of matrix
        :return: number of rows, number of columns
        :rtype: tuple
        """
        return self.mrow(), self.ncol()

    def isConformableAddSub(self, other):
        """
        Compares row and column lengths of two matrices to see if they
        are able to be added together
        :param other: second matrix to be added
        :type other: Matrix
        :return: T/F
        :rtype: bool
        """
        if self.m == other.m and self.n == other.n:
            return True
        return False

    def isConformableMult(self, other):
        """
        Compares row and column lengths of two matrices to see if they
        are able to be multiplied together via matrix multiplication
        :param other: second matrix to be added
        :type other: Matrix
        :return: T/F
        :rtype: bool
        """
        if self.n == other.m:
            return True
        return False

    def isSquareMat(self):
        """
        Are number of rows == number of columns?
        :return: T/F
        :rtype: bool
        """
        if self.m == self.n:
            return True
        return False

    def __add__(self, B):
        """
        Adds two matrices together, term by term,
        :param B: matrix to be added to self
        :type B: Matrix
        :return ret: summed matrix of self and B
        :rtype ret: Matrix
        """
        if self.isConformableAddSub(B):
            ret = Matrix(self.m, self.n)
            for x in range(self.mrow()):
                row = [sum(item) for item in zip(self[x], B[x])]
                ret[x] = row

            return ret
        else:
            raise MatrixError('Unable to add these two matrices')

    def __sub__(self, B):
        """Takes B (2nd matrix) as arg (Matrix object). Checks to see if self and B are conformable,
        then subtracts B from self if they are. Returns a new Matrix object with results.
        Error raised if self and B are not conformable."""
        if self.isConformableAddSub(B):
            ret = Matrix(self.m, self.n)
            for x in range(self.m):
                row = [item[0] - item[1] for item in zip(self.m[x], B[x])]
                ret[x] = row
            return ret
        else:
            raise MatrixError('Unable to add these two matrices')

    def __mul__(self, other):
        """Multiplies self to a scalar or another matrix.
        Takes other (Matrix object) as arg. Returns a new Matrix object with results."""
        # print(type(other))
        if isinstance(self, Matrix) and (isinstance(other, int) or isinstance(other, float)):
            return self.ScalarMultiply(other)
        elif isinstance(self, Matrix) and isinstance(other, Matrix):
            return self.MatrixMultiply(other)

    def ScalarMultiply(self, const):
        """Takes a const (int or float) as arg. Multiplies each element in Matrix by const.
        Returns a new Matrix object."""
        ret = Matrix(self.m, self.n)
        if isinstance(const, int) or isinstance(const, float):
            for x in range(self.m):
                ret.rows[x] = [item * const for item in self[x]]

            return ret
        else:
            raise MatrixError('Unable to multiply Matrix object with const')

    def MatrixMultiply(self, other):
        """Takes other (Matrix object) as arg. Checks to see if self and other are conformable.
        Multiplies them if they are. Returns a new Matrix object."""
        if self.isConformableMult(other):
            otherm, othern = other.getSize()
            other_t = other.getTranspose()
            mulmat = Matrix(self.m, othern)

            for x in range(self.m):
                for y in range(other_t.m):
                    mulmat[x][y] = sum([item[0] * item[1] for item in zip(self.rows[x], other_t[y])])

            return mulmat

        else:
            raise MatrixError('Matrices cannot be multiplied; self columns must equal other rows')

    def getTranspose(self):
        """Takes no args. Returns new, transposed Matrix object."""
        m, n = self.n, self.m
        mat = Matrix(m, n)
        mat.rows = [list(item) for item in zip(*self.rows)]

        return mat

    def interchangeRows(self, row1, row2):
        """Takes row1 (int) and row2 (int) as args, both of which designate a row in the Matrix.
        Switches the rows in place. Returns no values."""
        temp = self.rows[row1]
        self.rows[row1] = self.rows[row2]
        self.rows[row2] = temp

    def multiplyRowByConst(self, row, const):
        """Takes row (int) designated to by multiplied, and const(int or float) as args.
        Edits the row in place. Returns no values."""
        self.rows[row] *= const

    def addRows(self, row1, row2edit):
        """Takes row1 (int) and row2edit (int) as args, both of which designate a row in the Matrix.
        Edits row2edit in place. Returns no values."""
        r1, r2 = self.rows[row1], self.rows[row2edit]
        for i in r1:
            r1[i] += r2[i]

    def GaussJordanElimination(self, b):
        """This function reduces the augmented coeffecient matrix [A,b].
        Takes b (x by 1 matrix of coeffecients) as arg.
        # Unfinished."""
        E = 1
        C = Matrix(m=self.n, n=self.n + 1)

        list_of_col_val = []
        for j in range(C.m):
            list_of_col_val.append(C.m[j][0])
            print(list_of_col_val)

        p = max(list_of_col_val)

        for j in range(C.m):
            # compute pivot idx p, j<p<n exists in Cpj=max(Cij), j<=i<=n
            p = max([item for item in C[j]])
            print(p)


            if C.rows[p][j] == 0:
                E = 0
                break
            if p > j:
                C.interchangeRows(C.rows[p], C.rows[j])
            C.multiplyRowByConst(C.rows[j], 1 / C[j][j])
            # for i != j, - Cij * rj from ri

    def GaussElimination(self):
        """Takes no args. Supposed to return a matrix of rows, with the lower triangle as 0s.
        # Unfinished."""
        for j in range(min(self.n, self.m)):
            # pivot
            if self.rows[j][j] == 0:
                col = [self.rows[k][j] for k in range(j, self.m)]
                ipiv = col.index(max(col))
                temp = self.rows[j]
                self.rows[j] = self.rows[ipiv]
                self.rows[ipiv] = temp

            for i in range(j+1, self.m):
                ratio = self.rows[i][j] / self.rows[j][j]
                self.rows[i] = [self.rows[i][k] - (self.rows[j][k] * ratio) for k in range(self.n)]

        return self.rows


    def populateAugMat(self, A, b):
        """A and b are matrices augmented to give [A,b]. Takes A and b as args, returns new Matrix object.
        # Unfinished."""
        C = Matrix(m=self.n, n=self.n + 1)
        pass

    def getInverse(self):
        """Takes no args. Checks to see if square, then finds the inverse matrix.
        Also calls determinate finder. Returns new matrix object.
        # Unfinished."""
        pass

    def getDet(self, GR = True):
        """Takes GR (Gauss reduction method) as arg. If True, it follows the Gauss reduction algorithm.
        If False, method multiplies diagonals to get determinant.
        Finds the determinate of a square matrix.
        Returns determinant value (number- int or float)."""

        if GR is True:
            # Default. Follows Gauss reduction algorithm
            # Put Gauss Reduction form here
            print('not yet implemented')

            GE = self.GaussElimination()

            print(GE)

        elif GR is False:
            # Pseudo-recursive row and column eliminator
            if self.isSquareMat():
                if self.getSize() == (2, 2):
                    return (self.rows[0][0] * self.rows[1][1]) - (self.rows[0][1] * self.rows[1][0])
                else:
                    sum = 0
                    for x in range(len(self.rows[0])):
                        new = self.removeRowsColumns(i=self.rows[0], j=x).getDet(GR=False)
                        sum += (-1 ** (x % 2)) * self.rows[0][x] * new

                return sum
            else:
                raise MatrixError('Matrix must be a square matrix (n x n)')

    def removeRowsColumns(self, i, j):
        """Takes i and j (int) as args. i is the row, j is the column of the parent matrix whose
        indices will be removed. Returns newMat (Matrix object) of size m-1, n-1 with chosen
        row and column omitted"""
        newMat = Matrix(m=self.m - 1, n=self.n - 1)
        for x in range(len(self.rows)):
            for y in range(len(self.rows[x])):
                if x != i and y != j:
                    newMat.rows[x][y] = self.rows[x][y]
        return newMat

    def getAvgVec(self):
        """Takes no args. Finds average vector size by finding sum vector, then dividing by the number
        of rows used. Utilizes getSumVec(). Returns a 2x1 vector (Matrix object)."""
        summ = self.getSumVec()
        return summ * (1/self.m)

    def getSumVec(self):
        """Takes no args. Finds sum of rows of self (Matrix object).
        Returns a 2x1 vector (Matrix object)."""
        s1, s2 = 0, 0
        for i in self.rows:
            s1 += i[0]
            s2 += i[1]

        ret = Matrix(2, 1)
        ret[0][0] = s1
        ret[1][0] = s2

        return ret

    def getCovarMat(self, meanVec):
        """Takes meanVec (2x1 Matrix object) as arg. meanVec is derived from getAvgVec() above.
        Computes the covariance matrix for a given matrix of mx2 parent matrix.
        Returns a nxn matrix, the covariance matrix (Matrix object)."""
        # Setup
        diffVec = Matrix(2, 1)
        diff_list = []
        transpose_list = []
        nxnprod_list = []

        # Subtracting the mean vector from each of the measurement vectors in that class
        for row in range(len(self.rows)):
            for idx in range(2):
                diffVec[idx][0] = self.rows[row][idx] - meanVec.rows[idx][0]

            # Add to lists, and get n x n matrix
            diff_list.append(diffVec)
            transpose_list.append(diffVec.getTranspose())

        for vec_idx in range(len(diff_list)):
            nxnprod_list.append(diff_list[vec_idx] * transpose_list[vec_idx])

        # Accumulate sum
        total = Matrix(nxnprod_list[0].m, nxnprod_list[0].n)
        for nxn in nxnprod_list:
            total += nxn
        # Multiply by 1/k (k is # of vectors)
        return total * (1 / self.m)

    def getTrace(self):
        if self.isSquareMat():
            sum = 0
            for i in range(len(self.rows)):
                sum += self.rows[i][i]
            return sum
        else:
            raise MatrixError('Matrix must be a square matrix (n x n)')

    def findEigenvalues(self):
        if self.m == 2 and self.n == 2:
            return Smath.quadEqnSolver(a=1, b=-1*self.getTrace(), c=self.getDet(GR=False))
        else:
            print('Matrices > 2x2 coming soon')

    def findEigenvectors(self, lmd):
        # Want to find Lambda*I - A (a matrix)
        LImA = Matrix(m=self.m, n=self.n)
        for m in range(len(self.rows)):
            for n in range(len(self.rows[m])):
                # - A
                LImA.rows[m][n] = -1 * self.rows[m][n]
                # + Lambda*I
                if m == n:
                    LImA.rows[m][n] += lmd
        # print(LImA)

        # Solution vector
        sln_mat = Matrix(m=self.n, n=1)
        sln_mat.rows[0][0] = 1
        sln_mat.rows[1][0] = LImA.rows[0][0] / LImA.rows[0][1]

        # print('Solution matrix:', sln_mat, '##')

        ret = self * sln_mat
        return ret


    def getMagnitude(self):
        """Only works for 2x1 vector matrices. First value is x-coord, second is y-coord
        in Cartesian plane. Finds the hypotenuse of a right triangle with side lengths x and y."""
        if self.m == 2 and self.n == 1:
            return ((self.rows[0][0] ** 2) + (self.rows[1][0] ** 2)) ** 0.5
        else:
            raise MatrixError('Matrix must be a 2 x 1 vector')

    def getUnitVector(self):
        return self * (1/self.getMagnitude())