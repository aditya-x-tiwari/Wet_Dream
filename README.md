# Wet_Dream
Wet Dream is a proof-of-concept atmospheric water generator (AWG) designed to harvest water from ambient air, providing a sustainable source of drinking water. The project includes:

- Thermodynamic modeling of the vapour-compression refrigeration cycle
- Python scripts for calculating pressures, temperatures, enthalpies, and entropies at all key points in the cycle
- Coefficient of performance (COP) and efficiency analysis
- Potential for expansion to real-world air–water extraction prototypes

This repository is intended for engineers, students, and hobbyists interested in renewable water generation and thermodynamic cycle simulation.

"""
Advanced AWG Model with NTU-based Evaporator
Includes:
- Weather fetch
- Psychrometrics
- NTU-effectiveness evaporator model
- Full R134a cycle
- Airflow Reynolds check
- Parametric sweep
"""

import os
import numpy as np
import pandas as pd
import requests
from CoolProp.CoolProp import PropsSI
from dotenv import load_dotenv
import psychrolib

psychrolib.SetUnitSystem(psychrolib.SI)

# ---------------- USER INPUTS ----------------
load_dotenv()
OWM_API_KEY = os.getenv("OWM_API_KEY")
LOCATION = ("22.5726", "88.3639")  # Kolkata

refrigerant = "R134a"
EVAP_DEW_OFFSET = -10.0
COND_AMB_OFFSET = 3.0

# Heat exchanger assumptions
U = 80.0       # W/m2-K
A = 1.0        # m2
h_fg = 2.45e6  # J/kg

# Air properties
mu_air = 1.8e-5
D_tube = 0.01
flow_area = 0.5

m_vals = np.arange(0.1, 5.1, 0.1)
# ------------------------------------------------


def fetch_weather():
    if not OWM_API_KEY:
        raise ValueError("No API key found.")
    base = "https://api.openweathermap.org/data/2.5/weather"
    params = {"lat": LOCATION[0], "lon": LOCATION[1],
              "appid": OWM_API_KEY, "units": "metric"}
    r = requests.get(base, params=params, timeout=10)
    data = r.json()
    return (float(data["main"]["temp"]),
            float(data["main"]["pressure"]) * 100,
            float(data["main"]["humidity"]) / 100)


def isentropic_eff(PR):
    return np.clip(0.78 - 0.03 * max(0, PR - 2), 0.6, 0.85)


def run():

    T_amb, P_amb, RH = fetch_weather()
    T_dew = psychrolib.GetTDewPointFromRelHum(T_amb, RH)
    T_evap_air = T_dew + EVAP_DEW_OFFSET
    T_cond = T_amb + COND_AMB_OFFSET

    omega_in = psychrolib.GetHumRatioFromRelHum(T_amb, RH, P_amb)
    omega_out = psychrolib.GetHumRatioFromRelHum(T_evap_air, 1.0, P_amb)

    cp_moist = 1005 + omega_in * 1860
    rho_air = psychrolib.GetMoistAirDensity(T_amb, omega_in, P_amb)

    # Refrigerant side states
    T_evap_K = T_evap_air + 273.15
    T_cond_K = T_cond + 273.15

    P_evap = PropsSI("P", "T", T_evap_K, "Q", 1, refrigerant)
    P_cond = PropsSI("P", "T", T_cond_K, "Q", 0, refrigerant)

    h1 = PropsSI("H", "T", T_evap_K, "Q", 1, refrigerant)
    s1 = PropsSI("S", "T", T_evap_K, "Q", 1, refrigerant)
    h3 = PropsSI("H", "T", T_cond_K, "Q", 0, refrigerant)
    h4 = h3

    PR = P_cond / P_evap
    eta = isentropic_eff(PR)
    h2s = PropsSI("H", "P", P_cond, "S", s1, refrigerant)
    h2 = h1 + (h2s - h1) / eta

    q_evap = h1 - h4

    results = []

    for m_da in m_vals:

        C_air = m_da * cp_moist
        NTU = U * A / C_air if C_air > 0 else 0
        eps = 1 - np.exp(-NTU)

        # Required cooling
        T_target = T_dew - 5
        Q_s = m_da * cp_moist * (T_amb - T_target)
        water = max(0, m_da * (omega_in - omega_out))
        Q_l = water * h_fg
        Q_req = Q_s + Q_l

        # HX limited cooling
        Q_max = C_air * (T_amb - T_evap_air)
        Q_HX = eps * Q_max
        Q_actual = min(Q_req, Q_HX)

        # Refrigerant mass flow
        m_ref = Q_actual / q_evap if q_evap > 0 else 0
        W_comp = m_ref * (h2 - h1)
        COP = Q_actual / W_comp if W_comp > 0 else np.nan

        # Water production
        water_hr = water * 3600
        water_per_kWh = water_hr / (W_comp / 1000) if W_comp > 0 else np.nan

        # Reynolds number
        velocity = m_da / (rho_air * flow_area)
        Re = rho_air * velocity * D_tube / mu_air
        regime = "Laminar" if Re < 2300 else "Transitional" if Re < 4000 else "Turbulent"

        results.append([m_da, Q_actual, m_ref, W_comp,
                        COP, water_hr, water_per_kWh,
                        NTU, eps, Re, regime])

    cols = ["m_da", "Q_W", "m_ref", "W_comp",
            "COP", "water_kg_hr", "water_kg_per_kWh",
            "NTU", "effectiveness", "Re", "flow_regime"]

    df = pd.DataFrame(results, columns=cols)
    df.to_csv("awg_results.csv", index=False)

    best = df.loc[df["water_kg_per_kWh"].idxmax()]
    print("\nBest operating point:")
    print(best)

    return df


if __name__ == "__main__":
    df = run()
    print(df.head())