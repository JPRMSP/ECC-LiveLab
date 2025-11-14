# app.py
import streamlit as st
import hashlib
import random
import io
from math import isqrt
import matplotlib.pyplot as plt

st.set_page_config(layout="wide", page_title="ECC LiveLab — Extended")

# -----------------------
# Finite field helpers
# -----------------------
def inv_mod(k, p):
    # modular inverse via pow (p prime assumed)
    k = k % p
    if k == 0:
        raise ZeroDivisionError("Inverse of 0")
    return pow(k, p - 2, p)

def legendre_symbol(a, p):
    return pow(a, (p - 1) // 2, p)

def mod_sqrt(a, p):
    # Tonelli-Shanks would be general; for small primes we attempt brute force
    a %= p
    if a == 0:
        return 0
    if p % 4 == 3:
        x = pow(a, (p + 1) // 4, p)
        if (x * x) % p == a:
            return x
        return None
    # fallback brute force (p small for visualization/teaching)
    for x in range(p):
        if (x * x) % p == a:
            return x
    return None

# -----------------------
# Elliptic curve operations
# -----------------------
def is_on_curve(P, a, b, p):
    if P is None:
        return True
    x, y = P
    return (y * y - (x * x * x + a * x + b)) % p == 0

def ec_add(P, Q, a, p):
    if P is None: return Q
    if Q is None: return P
    x1, y1 = P
    x2, y2 = Q
    if x1 == x2 and (y1 + y2) % p == 0:
        return None
    if P == Q:
        if y1 == 0:
            return None
        m = (3 * x1 * x1 + a) * inv_mod(2 * y1, p) % p
    else:
        if (x2 - x1) % p == 0:
            return None
        m = (y2 - y1) * inv_mod(x2 - x1, p) % p
    x3 = (m * m - x1 - x2) % p
    y3 = (m * (x1 - x3) - y1) % p
    return (x3, y3)

def ec_mul(P, n, a, p):
    R = None
    Q = P
    while n > 0:
        if n & 1:
            R = ec_add(R, Q, a, p)
        Q = ec_add(Q, Q, a, p)
        n >>= 1
    return R

# -----------------------
# Curve enumeration & plotting
# -----------------------
def list_points(a, b, p):
    pts = []
    for x in range(p):
        rhs = (x*x*x + a*x + b) % p
        y = mod_sqrt(rhs, p)
        if y is not None:
            pts.append((x, y))
            if y != (p - y) % p:
                pts.append((x, p - y))
    return pts

def plot_curve(points, selected=None):
    xs = [pt[0] for pt in points]
    ys = [pt[1] for pt in points]
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.scatter(xs, ys, s=10)
    if selected:
        for label, pt in selected.items():
            if pt is not None:
                ax.scatter([pt[0]], [pt[1]], s=60, marker='x')
                ax.annotate(label, xy=pt, textcoords="offset points", xytext=(5,5))
    ax.set_xlabel("x")
    ax.set_ylabel("y")
    ax.set_title("Points on curve (toy visualization)")
    plt.tight_layout()
    return fig

# -----------------------
# Message <-> Point encoding (toy)
# - Map a short text to an integer then try to find x such that RHS is quadratic residue
# -----------------------
def msg_to_point(msg, a, b, p, max_tries=1000):
    h = int(hashlib.sha256(msg.encode()).hexdigest(), 16)
    for k in range(max_tries):
        x = (h + k) % p
        rhs = (x*x*x + a*x + b) % p
        y = mod_sqrt(rhs, p)
        if y is not None:
            return (x, y)
    raise ValueError("Failed to encode message to a point (try different curve or shorter message)")

def point_to_msg_dummy(P):
    # reverse mapping is not unique; we return a hex-like representation as "plaintext"
    if P is None:
        return ""
    return f"POINT({P[0]},{P[1]})"

# -----------------------
# EC-ElGamal
# -----------------------
def ec_elgamal_encrypt(M_point, G, Q, a, p):
    k = random.randrange(1, p-1)
    C1 = ec_mul(G, k, a, p)
    S = ec_mul(Q, k, a, p)  # kQ = k*d*G
    # "Add" message point and S
    C2 = ec_add(M_point, S, a, p)
    return (C1, C2), k

def ec_elgamal_decrypt(C, d, a, p):
    C1, C2 = C
    S = ec_mul(C1, d, a, p)
    # M = C2 - S -> add C2 + (-S)
    if S is None:
        return C2
    S_neg = (S[0], (-S[1]) % p)
    M = ec_add(C2, S_neg, a, p)
    return M

# -----------------------
# Toy "pairing" demo
# -----------------------
def toy_pairing_map(P, Q, p):
    # NOT a real pairing. This is an illustrative mapping that outputs a small integer derived
    # from both points so students can see that mapping group elements to another group
    # could leak relationships. DO NOT use as real pairing.
    if P is None or Q is None:
        return 0
    # simple deterministic map: combine x-coordinates and reduce mod p
    return (P[0] * 31 + Q[0] * 17 + P[1] * 7 + Q[1] * 13) % p

def pairing_breaks_ddh_demo(G, a_priv, b_priv, c_pub, G_point, a, p):
    # shows that e(aG, bG) ?= e(G, cG) for DDH decision via toy map
    A = ec_mul(G_point, a_priv, a, p)
    B = ec_mul(G_point, b_priv, a, p)
    C = ec_mul(G_point, c_pub, a, p)
    left = toy_pairing_map(A, B, p)
    right = toy_pairing_map(G_point, C, p)
    return left, right, left == right

# -----------------------
# Streamlit UI
# -----------------------
st.title("🔐 ECC LiveLab — Extended Playground")

col1, col2 = st.columns([2, 1])

with col1:
    st.header("Curve setup & visualization")
    p = st.number_input("Prime field p (small prime for demo & plotting)", value=97, min_value=3)
    a = st.number_input("Curve parameter a", value=2)
    b = st.number_input("Curve parameter b", value=3)
    # quick discriminant check
    disc = (4 * (a**3) + 27 * (b**2)) % p
    if disc == 0:
        st.error("Curve discriminant = 0 (singular). Pick different a/b.")
        st.stop()
    points = list_points(a, b, p)
    st.write(f"Number of affine points found (including ±y): {len(points)} (toy)")
    st.pyplot(plot_curve(points))

    st.subheader("Choose generator point G")
    Gx = st.number_input("Gx", value=3)
    Gy = st.number_input("Gy", value=6)
    G = (Gx, Gy)
    if not is_on_curve(G, a, b, p):
        st.warning("G is not on the curve (choose another G).")
    else:
        st.success("G is on the curve (toy check).")

    st.subheader("Scalar multiplication demo")
    k = st.number_input("Scalar k", value=7)
    kG = ec_mul(G, k, a, p) if is_on_curve(G, a, b, p) else None
    st.write("k * G =", kG)

with col2:
    st.header("Key generation & ECDH")
    d = st.number_input("Private key d (choose small integer)", value=11)
    Q = ec_mul(G, d, a, p) if is_on_curve(G, a, b, p) else None
    st.write("Public key Q =", Q)
    if st.button("Download keypair"):
        priv_txt = f"PRIVATE KEY (toy)\ncurve: y^2 = x^3 + {a}x + {b} over F_{p}\nprivate: {d}\n"
        pub_txt = f"PUBLIC KEY (toy)\ncurve: y^2 = x^3 + {a}x + {b} over F_{p}\npublic: {Q}\n"
        st.download_button("Download PRIVATE", priv_txt, file_name="ecc_private.txt")
        st.download_button("Download PUBLIC", pub_txt, file_name="ecc_public.txt")

st.markdown("---")

colA, colB = st.columns(2)

with colA:
    st.header("ECDH Demo (Alice & Bob)")
    a_secret = st.number_input("Alice secret a", value=5, key="a_secret")
    b_secret = st.number_input("Bob secret b", value=9, key="b_secret")
    A_pub = ec_mul(G, a_secret, a, p)
    B_pub = ec_mul(G, b_secret, a, p)
    shared1 = ec_mul(B_pub, a_secret, a, p)
    shared2 = ec_mul(A_pub, b_secret, a, p)
    st.write("Alice public A =", A_pub)
    st.write("Bob public B =", B_pub)
    st.write("Shared (Alice) =", shared1)
    st.write("Shared (Bob)   =", shared2)
    if shared1 == shared2:
        st.success("Shared keys match (ECDH OK).")
    else:
        st.error("Shared keys do not match (check inputs).")

with colB:
    st.header("ECDSA (toy) — quick sign/verify")
    msg = st.text_input("Message to sign (short)", "hello ecc")
    n = p  # toy group order substitute
    if st.button("Sign message"):
        # toy signing similar to earlier example
        h = int(hashlib.sha256(msg.encode()).hexdigest(), 16)
        k_rand = (h + d) % n or 1
        R = ec_mul(G, k_rand, a, p)
        if R is None:
            st.error("Bad ephemeral k produced R=∞; try different params")
        else:
            r = R[0] % n
            s = (inv_mod(k_rand, n) * (h + r * d)) % n
            st.write("Signature (r, s):", (r, s))
            # verify
            w = inv_mod(s, n)
            u1 = (h * w) % n
            u2 = (r * w) % n
            P = ec_add(ec_mul(G, u1, a, p), ec_mul(Q, u2, a, p), a, p)
            valid = (P is not None) and (P[0] % n == r)
            st.write("Verification result:", valid)

st.markdown("---")

st.header("EC-ElGamal Encryption (toy)")

colE1, colE2 = st.columns(2)
with colE1:
    plaintext = st.text_input("Plaintext (short)", "HELLO")
    try:
        M_pt = msg_to_point(plaintext, a, b, p)
        st.write("Encoded message point M =", M_pt)
    except Exception as e:
        st.error("Encoding failed: " + str(e))
        M_pt = None

with colE2:
    if st.button("Encrypt with recipient public Q"):
        if M_pt is None:
            st.error("No valid message point")
        elif Q is None:
            st.error("No recipient public key Q")
        else:
            (C1, C2), k_used = ec_elgamal_encrypt(M_pt, G, Q, a, p)
            st.write("Ciphertext C1, C2 =", (C1, C2))
            st.write("Ephemeral k used (toy):", k_used)
            # offer download
            ct_txt = f"C1: {C1}\nC2: {C2}\n"
            st.download_button("Download ciphertext", ct_txt, file_name="ciphertext.txt")

if st.button("Decrypt last ciphertext (toy)"):
    try:
        C1, C2
    except NameError:
        st.error("No ciphertext present in this session; encrypt first")
    else:
        M_rec = ec_elgamal_decrypt((C1, C2), d, a, p)
        st.write("Recovered point M' =", M_rec)
        st.write("Recovered data (toy repr) =", point_to_msg_dummy(M_rec))

st.markdown("---")

st.header("Toy Pairing Demo (illustrative only)")

st.info("This is a simple mapping (NOT a real Weil/Tate pairing). It's a classroom demo to show how mapping points to another group can leak relationships and thus break DDH on small toy curves.")

c_secret = st.number_input("Third secret c (for DDH challenge)", value=7, key="c_secret")
left, right, eq = pairing_breaks_ddh_demo(G, a_secret, b_secret, c_secret, G, a, p)
st.write("e(aG, bG)  (toy map) =", left)
st.write("e(G, cG)   (toy map) =", right)
if eq:
    st.success("Mapping values equal — toy-demo suggests DDH instance is 'yes'")
else:
    st.warning("Mapping values differ — toy-demo suggests DDH instance is 'no'")

st.markdown("---")
