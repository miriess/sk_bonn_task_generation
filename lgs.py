"""Lineare Algebra Aufgaben
"""


from sympy import *


init_session()


def makeLGS(dim, entryrange=10, resultrange=5, maxdet=50):
    """Generiert in dim x dim LGS mit Einträgen zwischen -entryrange und
    +entryrange sowie Ergebnissen zwischen -resultrange und +resultrange.
    """

    from random import randint

    A = Matrix([[randint(-entryrange, entryrange) for i in range(dim)]
                for i in range(dim)])

    R = Matrix([randint(-resultrange, resultrange) for i in range(dim)])

    if A.det() == 0:
        return makeLGS(dim, entryrange, resultrange, maxdet)
    elif abs(A.det()) > maxdet:
        return makeLGS(dim, entryrange, resultrange, maxdet)
    else:
        return [A, A*R, R, A.det()]


def linAbh(clu, cla, dim, com=0, coordrange=9, paramrange=3):
    """Erzeugt eine Liste von Spaltenvektoren der Dimension dim, von denen clu
    linear unabhängig und cla linear abhängig sind -- Einträge im Betrag <=
    coordrange und die Linearkombinationsparameter sind <= paramrange.
    """

    from random import choice

    coordL = (
        [i+1 for i in range(coordrange)] + [0] +
        [-i-1 for i in range(coordrange)]
        )

    paramL = (
        [i+1 for i in range(paramrange)] + [0] +
        [-i-1 for i in range(paramrange)]
        )

    U = Matrix([choice(coordL) for i in range(dim)])

    for i in range(clu-1):
        U = U.col_insert(i+1, Matrix([choice(coordL) for j in range(dim)]))
        while U.rank() <= i+1:
            U.col_del(i+1)
            U = U.col_insert(i+1, Matrix([choice(coordL) for j in range(dim)]))

    for i in range(cla):
        v = Matrix([0 for j in range(dim)])
        if com:
            for j in range(U.shape[1]):
                v = v + choice(paramL) * U.col(j)
        else:
            for j in range(clu):
                v = v + choice(paramL) * U.col(j)
        U = U.col_insert(U.shape[1], v)

    result = [U.col(i) for i in range(U.shape[1])]

    return result
