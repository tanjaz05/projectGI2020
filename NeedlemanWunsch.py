import numpy

class NeedlemanWunsch:

    def __init__(self, match, mismatch, insertion):
        self.match = match
        self.mismatch = mismatch
        self.insertion = insertion

    def scoringMatrix(self, a, b):
        if a == b: return self.match
        if a == '_' or b == '_': return self.insertion
        return self.mismatch

    def globalAlignment(self, x, y, s):
        D = numpy.zeros((len(x) + 1, len(y) + 1), dtype=int)

        for i in range(1, len(x) + 1):
            D[i, 0] = D[i - 1, 0] + s(x[i - 1], '_')
        for j in range(1, len(y) + 1):
            D[0, j] = D[0, j - 1] + s('_', y[j - 1])

        for i in range(1, len(x) + 1):
            for j in range(1, len(y) + 1):
                D[i, j] = max(D[i - 1, j] + s(x[i - 1], '_'),
                              D[i, j - 1] + s('_', y[j - 1]),
                              D[i - 1, j - 1] + s(x[i - 1], y[j - 1]))
        return D, D[len(x), len(y)]

    def getTranscript(self, D, x, y):
        t = ''
        i = len(x)
        j = len(y)
        while i != 0 and j != 0:
            delta = 1 if x[i - 1] != y[j - 1] else 0
            Del = D[i - 1, j] + 1
            I = D[i, j - 1] + 1
            M = D[i - 1, j - 1] + delta

            # diagonal was the best
            if M <= I and M <= Del:
                # for diagonal we check wether characters on that position were the same of different
                # based on that we add M (match) or R (replacement) edit operation
                t += 'M' if delta == 0 else 'R'
                i = i - 1
                j = j - 1

            # horizontal was the best
            elif I <= Del:
                t += 'I'
                j = j - 1

            # vertical was the best
            else:
                t += 'D'
                i = i - 1

        if i != 0:
            t += 'D' * i

        if j != 0:
            t += 'I' * j

        # we revert string in order to get edit transcript
        return t[::-1]

    def traceback(self, x, y, V, s):
        # initializing starting position cell(n,m)
        i = len(x)
        j = len(y)

        # initializing strings we use to represent alignments in x, y, edit transcript and global alignment
        ax, ay, am, tr = '', '', '', ''

        # exit condition is when we reach cell (0,0)
        while i > 0 or j > 0:

            # calculating diagonal, horizontal and vertical scores for current cell
            d, v, h = -100, -100, -100

            if i > 0 and j > 0:
                delta = 1 if x[i - 1] == y[j - 1] else 0
                d = V[i - 1, j - 1] + s(x[i - 1], y[j - 1])  # diagonal movement
            if i > 0: v = V[i - 1, j] + s(x[i - 1], '_')  # vertical movement
            if j > 0: h = V[i, j - 1] + s('_', y[j - 1])  # horizontal movement

            # backtracing to next (previous) cell
            if d >= v and d >= h:
                ax += x[i - 1]
                ay += y[j - 1]
                if delta == 1:
                    tr += 'M'
                    am += '|'
                else:
                    tr += 'R'
                    am += ' '
                i -= 1
                j -= 1
            elif v >= h:
                ax += x[i - 1]
                ay += '_'
                tr += 'D'
                am += ' '
                i -= 1
            else:
                ay += y[j - 1]
                ax += '_'
                tr += 'I'
                am += ' '
                j -= 1

        alignment = '\n'.join([ax[::-1], am[::-1], ay[::-1]])
        return alignment, tr[::-1]
