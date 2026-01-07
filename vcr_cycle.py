import CoolProp.CoolProp as CP

# -----------------------------
# INPUT DATA
# -----------------------------
fluid = 'R134a'

Tevap = 5 + 273.15      # Evaporator temperature [K]
Tcond = 40 + 273.15     # Condenser temperature [K]

# -----------------------------
# STATE 1: Saturated vapor at evaporator
# -----------------------------
P1 = CP.PropsSI('P', 'T', Tevap, 'Q', 1, fluid)
T1 = Tevap
h1 = CP.PropsSI('H', 'P', P1, 'Q', 1, fluid)
s1 = CP.PropsSI('S', 'P', P1, 'Q', 1, fluid)

# -----------------------------
# STATE 2: Isentropic compression
# -----------------------------
P2 = CP.PropsSI('P', 'T', Tcond, 'Q', 0, fluid)
s2 = s1
h2 = CP.PropsSI('H', 'P', P2, 'S', s2, fluid)
T2 = CP.PropsSI('T', 'P', P2, 'S', s2, fluid)

# -----------------------------
# STATE 3: Saturated liquid at condenser
# -----------------------------
P3 = P2
T3 = Tcond
h3 = CP.PropsSI('H', 'P', P3, 'Q', 0, fluid)
s3 = CP.PropsSI('S', 'P', P3, 'Q', 0, fluid)

# -----------------------------
# STATE 4: Throttling (h4 = h3)
# -----------------------------
P4 = P1
h4 = h3
T4 = CP.PropsSI('T', 'P', P4, 'H', h4, fluid)
s4 = CP.PropsSI('S', 'P', P4, 'H', h4, fluid)

# -----------------------------
# PERFORMANCE CALCULATIONS
# -----------------------------
q_evap = h1 - h4
w_comp = h2 - h1
COP = q_evap / w_comp

# -----------------------------
# PRINT RESULTS
# -----------------------------
print("\nVAPOUR COMPRESSION CYCLE RESULTS\n")

print("State 1 (Evaporator outlet):")
print(f"P = {P1/1e5:.2f} bar, T = {T1-273.15:.2f} °C, h = {h1/1000:.2f} kJ/kg, s = {s1/1000:.4f} kJ/kg.K\n")

print("State 2 (Compressor outlet):")
print(f"P = {P2/1e5:.2f} bar, T = {T2-273.15:.2f} °C, h = {h2/1000:.2f} kJ/kg, s = {s2/1000:.4f} kJ/kg.K\n")

print("State 3 (Condenser outlet):")
print(f"P = {P3/1e5:.2f} bar, T = {T3-273.15:.2f} °C, h = {h3/1000:.2f} kJ/kg, s = {s3/1000:.4f} kJ/kg.K\n")

print("State 4 (Expansion valve outlet):")
print(f"P = {P4/1e5:.2f} bar, T = {T4-273.15:.2f} °C, h = {h4/1000:.2f} kJ/kg, s = {s4/1000:.4f} kJ/kg.K\n")

print(f"Refrigeration Effect = {q_evap/1000:.2f} kJ/kg")
print(f"Compressor Work      = {w_comp/1000:.2f} kJ/kg")
print(f"COP of Cycle         = {COP:.2f}")
