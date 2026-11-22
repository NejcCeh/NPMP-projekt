import numpy as np

# IZRAČUN ELEKTROSTATIČNE ENERGIJE NOTRANJE CELICE M PRI RAZLIČNIH GEOMETRIJSKIH ODMIKIH V QCA VEČINSKIH VRATIH

# fizikalne konstante
e = 1.602176634e-19 # velikost osnovnega električnega naboja
eps_0 = 8.8541878128e-12 # dielektričnost vakuuma
eps_r = 12.9 # relativna dielektričnost materiala

# geometrija QCA celic [nm]
a = 9 # razdalja med središči sosednjih kvantnih pik znotraj celice
d0 = 20 # razdalja med središčema sosednjih QCA celic v osnovni konfiguraciji

# kvantne pike
dots = np.array([
    [-a/2,  a/2], # 1: zgoraj levo
    [ a/2,  a/2], # 2: zgoraj desno
    [ a/2, -a/2], # 3: spodaj desno
    [-a/2, -a/2]  # 4: spodaj levo
])

# položaji vhodnih celic [nm]
centers = {
    "A": np.array([0, d0]),
    "B": np.array([-d0, 0]),
    "C": np.array([0, -d0])
}

# izbrana vhodna kombinacija
polarizations = {
    "A": 1,
    "B": -1,
    "C": -1
}


def occupied_dots(P):
    """Vrne indekse pik, v katerih sta elektrona."""
    if P == 1:
        return [0, 2] # indeksi prve in tretje pike
    return [1, 3]


def Coulomb_energy(r_nm):
    """Izračuna elektrostatično energijo med dvema elektronoma [J]."""
    r_m = r_nm * 1e-9 # r_nm...razdalja med elektronoma v nm
    return e**2 / (4 * np.pi * eps_0 * eps_r * r_m)


def total_energy(displacement, P_M):
    """Sprejme odmik notranje celice in njeno polarizacijo ter vrne skupno interakcijsko energijo vhodov z M."""
    center_M = np.array([0.0, displacement])
    energy = 0

    for name, center in centers.items(): # zanka po vhodnih celicah
        P_input = polarizations[name]

        for i in occupied_dots(P_input): # vsak od dveh elektronov v vhodni celici
            for j in occupied_dots(P_M): # za vsakega od njiju izbere oba elektrona notranje celice
                pos_input = center + dots[i] # abs. koord. elektrona vhodne celice
                pos_M = center_M + dots[j] # abs. koord. elektrona notranje celice

                r = np.linalg.norm(pos_input - pos_M) # razdalja med elektronoma
                energy += Coulomb_energy(r)

    # J -> meV
    return energy / e * 1000


print("Odmik | E(+1) | E(-1) | E(+1)-E(-1) | Napoved")
print("-" * 66)

for displacement in [0.0, 1.0, 2.0]:
    E_plus = total_energy(displacement, 1)
    E_minus = total_energy(displacement, -1)

    delta_E = E_plus - E_minus

    predicted_P = 1 if delta_E < 0 else -1

    print(
        f"{displacement:4.1f} nm | "
        f"{E_plus:7.3f} | "
        f"{E_minus:7.3f} | "
        f"{delta_E:+10.3f} | "
        f"{predicted_P:+d}"
    )