"""Stellt Funktionen für Analysisaufgaben zur Verfügung."""


from sympy import init_session
init_session()


def spInt():
    """Berechnet eine Schnittpunktaufgabe mit Parabel und Gerade.
    Inklusive Fläche!
    Berichtet eine Liste [p,g,sp,int] mit
    p Parabel
    g Gerade
    sp Liste der Schnittpunkte
    int Werte des Integrals
    """

    from random import randint
    from sympy import symbols, expand, integrate, S

    x = symbols('x')

    nmbs = [randint(-6, 6) / S(3), randint(-10, 10), randint(-10, 10),
            randint(-10, 10), randint(-10, 10)]

    p = expand(nmbs[0]*(x-nmbs[1])*(x-nmbs[2])+(nmbs[3]*x + nmbs[4]))
    g = nmbs[3]*x+nmbs[4]

    if abs(integrate(p-g, (x, nmbs[1], nmbs[2]))) < 200:
        return [p, g, [nmbs[1], nmbs[2]],
                abs(integrate(p-g, (x, nmbs[1], nmbs[2])))]
    else:
        return spInt()


def findcost(X1, X2, X3, A):
    """Finde eine Preisuntergrenze zu einer Gewinnfunktion, so dass:
    - Die Gewinnfunktion ganzzahlige Nullstellen aus den Listen
        X1, X2, X3 aufweist.
    - Die resultierende Kostenfunktion monoton steigend ist.
    - Der führende Koeffizient aus der Liste A kommt.
    """

    from sympy import symbols, solve

    x = symbols('x')

    TestVar = []

    for x1 in X1:
        for x2 in X2:
            for x3 in X3:
                for a in A:
                    TestVar.append((x1, x2, x3, a))

    erg = []

    for set in TestVar:
        print(set)
        G = -set[3]*(x-set[0])*(x-set[1])*(x-set[2])
        WP = solve(G.diff(x, 2), x)[0]
        mWP = G.diff(x).subs(x, WP)
        erg.append((G, WP, mWP))

    return erg


def vollstK(G, p, varstr='x'):
    """Berechnet die Kennzahlen für vollständige Konkurrenz und gibt sie
    als Liste aus.
    """

    from sympy import symbols, solve, im

    x = symbols(varstr)

    E = p * x
    K = (E-G).expand()

    print('Gewinnfunktion ist ' + str(G.expand()))
    print('Erlösfunktion ist ' + str(E))
    print('Kostenfunktion ist ' + str(K))
    print('Grenzkostenfunktion ist ' + str(K.diff(x).expand()))
    print('Nullstellen der Gewinnfunktion sind ' + str(solve(G, x)))
    print('...gerundet ist das ' + str([a.evalf(6) for a in solve(G, x)]))

    Gmax = max(solve(G.diff(x), x))

    print('Gewinnmaximum von {0:.3f} bei x = {1:.3f}'.format(
        G.subs(x, Gmax).evalf(), Gmax.evalf()))

    gKmin = max(solve(K.diff(x, 2), x))

    print('Grenzkostenminimum von {0:.3f} bei x = {1:.3f}'.format(
        K.diff(x).subs(x, gKmin).evalf(), gKmin.evalf()))

    kV = ((K - K.subs(x, 0))/x).expand()
    k = (K/x).expand()

    print('Variable Stückkostenfunktion ist ' + str(kV))
    print('Stückkostenfunktion ist ' + str(k))

    Bmin = max(solve(kV.diff(x), x))

    print('Betriebsminimum mit variablen Stückkosten von  {0:.3f} bei x = {1:.3f}'.format(
        kV.subs(x, Bmin).evalf(), Bmin.evalf()))

    Bopt = max([a for a in solve(k.diff(x), x) if im(a) == 0])

    print('Betriebsoptimum mit Stückkosten von  {0:.3f} bei x = {1:.3f}'.format(
        k.subs(x, Bopt).evalf(), Bopt.evalf()))


def getVK(Xmin, gKmin, p, gZ, Bopt):
    """Berechnet eine Aufgabe aus dem Bereich der vollst. Konkurrenz mit den
    den Daten Xmin = x-Wert des Grenzkostenminimums, gKmin = Steigung von K(x)
    in Xmin, p = Preis des Produkts, gZ einer Grenze der Gewinnzone
    und Bopt = Betriebsoptimum.
    """

    from sympy import solve, symbols

    a, b, c, d, x = symbols('a b c d x')

    erg = solve([6*a*Xmin + 2*b,
                3*a*Xmin**2 + 2 * b * Xmin + c - gKmin,
                p * gZ - (a*gZ**3 + b*gZ**2 + c*gZ + d),
                2*a*Bopt**3 + b*Bopt**2 - d], [a, b, c, d])

    K = erg[a] * x**3 + erg[b] * x**2 + erg[c] * x + erg[d]

    E = p * x

    G = (E - K).expand()

    vollstK(G, p)

    return [K, E, G]


def getPolyp(Xmin, gKmin, b0, m, gZ, Bopt, report=0):
    """Berechnet eine Aufgabe im Polypol."""

    from sympy import solve, symbols

    a, b, c, d, x = symbols('a b c d x')

    erg = solve([
        6*a*Xmin + 2*b,
        3*a*Xmin**2 + 2 * b * Xmin + c - gKmin,
        (m*gZ + b0) * gZ - (a*gZ**3 + b*gZ**2 + c*gZ + d),
        2*a*Bopt**3 + b*Bopt**2 - d
        ], [a, b, c, d])

    K = erg[a] * x**3 + erg[b] * x**2 + erg[c] * x + erg[d]

    p = m * x + b0

    E = (p * x).expand()

    G = (E - K).expand()

    if report:
        from analysis import analysePolyp
        analysePolyp([K, p, E, G])

    return [K, p, E, G]


def analysePolyp(result):

    from sympy import solve, symbols, im

    x = symbols('x')

    K = result[0]

    p = result[1]

    E = p * x

    G = E - K

    print('Gewinnfunktion ist ' + str(G.expand()))
    print('Erlösfunktion ist ' + str(E))
    print('Kostenfunktion ist ' + str(K))
    print('Grenzkostenfunktion ist ' + str(K.diff(x).expand()))
    print('Nullstellen der Gewinnfunktion sind ' + str(solve(G, x)))
    print('...gerundet ist das ' + str([a.evalf(6) for a in solve(G, x)]))

    Gmax = max(solve(G.diff(x), x))

    print('Gewinnmaximum von {0:.3f} bei x = {1:.3f}'.format(
        G.subs(x, Gmax).evalf(), Gmax.evalf()))

    gKmin = max(solve(K.diff(x, 2), x))

    print('Grenzkostenminimum von {0:.3f} bei x = {1:.3f}'.format(
        K.diff(x).subs(x, gKmin).evalf(), gKmin.evalf()))

    kV = ((K - K.subs(x, 0))/x).expand()
    k = (K/x).expand()

    print('Variable Stückkostenfunktion ist ' + str(kV))
    print('Stückkostenfunktion ist ' + str(k))

    Bmin = max(solve(kV.diff(x), x))

    print('Betriebsminimum mit variablen Stückkosten von  {0:.3f} bei x = {1:.3f}'.format(
        kV.subs(x, Bmin).evalf(), Bmin.evalf()))

    Bopt = max([a for a in solve(k.diff(x), x) if im(a) == 0])

    print('Betriebsoptimum mit Stückkosten von  {0:.3f} bei x = {1:.3f}'.format(
        k.subs(x, Bopt).evalf(), Bopt.evalf()))
