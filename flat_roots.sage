#sage10.3 Halois
# if ur running a brute force attack on a polynomial, this is not your weapon, Too many system calls can drain your time


DEBUG = False  # set on True to see warnings 



def flat_roots(f, bounds, m=1, d=None):
    import itertools, subprocess
    if DEBUG:
        import warnings
        warnings.filterwarnings('ignore')
    def dump(M):
        return "[{}]".format("\n".join("[{}]".format(" ".join(map(str, r))) for r in M))

    def parse_row(r):
        return map(ZZ, r.split())
    
    def parse(x):
        x = x.replace("\n", "")
        assert x[0] == "["
        assert x[-1] == "]"
        return Matrix(ZZ, map(parse_row, x[2:-2].split("][")))
        
    def flatter(M):
        try:
            PATH = FLATTER_PATH
        except NameError:
            print("Warning: using `flatter` if it can be found on $PATH. Set `FLATTER_PATH` to the command if this fails.") if DEBUG else None
            PATH = "flatter"

        try:
            ARGS = FLATTER_ARGS
        except:
            ARGS = []

        return parse(subprocess.check_output([PATH, *ARGS], input=dump(M).encode()).decode())


    # Because univariate and multivariate are just different enough to make you crazy , but in sage10.3 we get a grace way
    PR = PolynomialRing(f.parent(), 0, ()).flattening_morphism().codomain()
    f = PR(f)
    
    orig = f
    if not d:
        d = f.degree()

    R = f.base_ring()
    N = R.cardinality()
    
    f /= f.coefficients().pop(0)
    f = f.change_ring(ZZ)

    G = Sequence([], f.parent())
    for i in range(m+1):
        base = N^(m-i) * f^i
        for shifts in itertools.product(range(d), repeat=len(f.parent().gens())):
            g = base * prod(map(power, f.parent().gens(), shifts))
            G.append(g)

    B, monomials = G.coefficient_matrix()
    monomials = vector(monomials)

    factors = [monomial(*bounds) for monomial in monomials]
    for i, factor in enumerate(factors):
        B.rescale_col(i, factor)

    B = flatter(B.dense_matrix())

    B = B.change_ring(QQ)
    for i, factor in enumerate(factors):
        B.rescale_col(i, 1/factor)

    H = Sequence([], f.parent().change_ring(QQ))
    roots = []
    for h in filter(None, B*monomials):
        H.append(h)
        I = H.ideal()
        if I.dimension() == -1:
            H.pop()
        elif I.dimension() == 0:
            roots = []
            for root in I.variety(ring=ZZ):
                root = tuple(R(root[var]) for var in f.parent().gens())
                if gcd(ZZ(orig(*root)), N) != 1:
                    roots.append(root)
            return roots
    return roots


def install_flatter():
    import os
    if not os.path.exists("flatter/bin/flatter"):
        os.system("rm -rf flatter \
                  && git clone https://github.com/keeganryan/flatter \
                  && cd flatter && sed -i 's/SHARED/STATIC/g' src/CMakeLists.txt \
                  && mkdir build \
                  && cd build \
                  && cmake .. \
                    -DCMAKE_BUILD_TYPE=Release \
                    -DCMAKE_INSTALL_PREFIX=..\
                    -DCMAKE_CXX_FLAGS='-march=native' \
                  && make install")
    global FLATTER_PATH
    FLATTER_PATH = os.path.join(os.getcwd(), "flatter", "bin", "flatter")