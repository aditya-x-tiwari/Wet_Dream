"""

What it does:

## Adaptive to Different Weather Conditions 
- No realtime sensors required, just internet connection or can also be input manually too.
- Fetches ambient temperature, pressure, and relative humidity from OpenWeatherMap (or fallback to manual input).
- Computes psychrometric properties (dew point, humidity ratio, moist-air enthalpy).

## Thermodynamic Calculations
- For a standard vapour compression cycle assuming no subcooling or superheating the specific enthalpies are found which will be used later this is done via CoolProp.
-  Libraries like Psychrolib are used to find various characterstics of refrigerant and air at different humid conditions.

## Heat Transfer Calculations
- For the efficient water generation we must not go below the dew point theoretically, (we have taken a tolernace of (+5)) thus the optimised amount of heat for a given flow rate of air is determined (including sensible + latent).
- We also need to determine the flow characterstics of air hence Reynolds and Nusselt Numbers are determined to identify flow characterstics.
- Prior to going towards the refrigerant side we have computed the thermal resistance and area and thereby computing NTU and finally effectiveness so that we can scale down the possible heat flow.
- Now, For a range of refrigerant mass flow rates: 0.1 to 5.0 kg/s (per step), we compute computes:
    * theoretical refrigerant-side cooling Q_ref
    * compressor power and COP
    * Reynolds and Nusselt Numbers to determine flow characterstics
    * water production rate (kg/hr) and water per kWh

## Results    
- A CSV with all the data is printed and finally we have calculated the the most optimum refrigerant flow required.
- Finally all of the data is made to be visualised for laymans and thus the program comes to end doing its job. 


Assumptions (editable in code):
- Evaporator refrigerant saturation temperature = (dew_point - 10°C)
- Condenser refrigerant saturation temperature = (T_ambient + 3°C)
- Evaporator treated as a heat exchanger with refrigerant phase-change side (C_r ~= 0)
- Isentropic compressor efficiency is estimated by a simple function (editable)
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
OWM_API_KEY = os.getenv("OWM_API_KEY") # not exposing api key directly as it's confidential.
print(OWM_API_KEY)
if not OWM_API_KEY:
    print("⚠️ No OpenWeatherMap API key found in environment variables.")

# Provide location: either (lat, lon) tuple or city name string
# Examples:
# LOCATION = ("22.5726","88.3639")   # lat,lon for Kolkata
# LOCATION = "Kolkata,IN"    
LOCATION = ("22.5726", "88.3639")  # Kolkata

# -------- PHYSICAL / DESIGN DEFAULTS (editable) ----------

# Refrigerant properties
refrigerant = "R134a"
# Temperatures offsets (°C) these are usually determined based on the model
EVAP_DEW_OFFSET = -10.0     # T_evap_air = dew_point + EVAP_DEW_OFFSET (note: offset is negative here)
COND_AMB_OFFSET = +3.0      # T_cond = T_ambient + COND_AMB_OFFSET

# ---- Sweep refrigerant mass flow ----
m_ref_vals = np.round(np.arange(0.1, 20.1, 0.1),1)


# Heat exchanger assumptions
U = 80.0       # W/m2-K # this contains a lot of other variables under consideration too but we have simplified things up for the time being.
A = 1.0        # m2 # this also involve complexities like it will be cross sectional area of fin that is circle with r as radius * no.of evaporator tubes in each fin *2 (two sides of fin * number of total fins in one evaporator * number of evaporaters
L_vap = 2.45e6  # J/kg latent heat of vapourization of water its a constant
D_tube = 0.01
N_evaporator= 1 # number of evaporator in parallel
A_tube=(np.pi * D_tube**2 / 4) * N_evaporator   # total internal flow area




# ------------------------------------------------


def fetch_weather():
    """
    Tries OpenWeatherMap.
    If no API key or failure → switches to manual input.
    """

    # ---- Try API first ----
    if OWM_API_KEY:
        try:
            base = "https://api.openweathermap.org/data/2.5/weather"
            params = {
                "lat": LOCATION[0],
                "lon": LOCATION[1],
                "appid": OWM_API_KEY,
                "units": "metric"
            }
            r = requests.get(base, params=params, timeout=10)
            r.raise_for_status()
            data = r.json()

            print("Weather fetched from OpenWeatherMap.")

            return (
                float(data["main"]["temp"]),
                float(data["main"]["pressure"]) * 100,
                float(data["main"]["humidity"]) / 100
            )

        except Exception as e:
            print("API fetch failed:", e)

    # ---- Manual fallback ----
    print("\nUsing manual weather input.")

    try:
        filename = 'input_data.csv'

        with open(filename, mode='r') as file:
            reader = csv.DictReader(file)

            # Skip header row (not needed for DictReader actually)
            # next(reader)

            for row in reader:
                T = float(row['Ambient Temperature (°C)'])
                P = float(row['Ambient Pressure (Pa)'])
                RH = float(row['Relative Humidity (0–1)'])
                return T, P, RH

    except Exception:
        print("Input failed — using safe default values.")
        # Default: Kolkata typical summer
        return 32.0, 101325.0, 0.7
def isentropic_eff(PR):
    return np.clip(0.78 - 0.03 * max(0, PR - 2), 0.6, 0.85)


def run():

    T_amb, P_amb, RH = fetch_weather()
    T_dew = psychrolib.GetTDewPointFromRelHum(T_amb, RH)
    T_evap_air = T_dew + EVAP_DEW_OFFSET
    T_cond = T_amb + COND_AMB_OFFSET


    print(f"Ambient: {T_amb:.2f} °C, RH: {RH*100:.1f} %, Dew point: {T_dew:.2f} °C")
    print(f"Assumed evaporator air temp (out): {T_evap_air:.2f} °C")
    print(f"Assumed condenser saturation temp: {T_cond:.2f} °C")



    # Precompute enthalpies for moist air inlet/outlet at a given humidity ratio
    omega_in = psychrolib.GetHumRatioFromRelHum(T_amb, RH, P_amb)  # kg water/kg of inlet air
    omega_out = psychrolib.GetHumRatioFromRelHum(T_evap_air, 1.0, P_amb)  # kg water/kg of outlet air 

    h_in = psychrolib.GetMoistAirEnthalpy(T_amb, omega_in)   # J/kg dry air
    h_out_saturated = psychrolib.GetMoistAirEnthalpy(T_evap_air, omega_out)

    cp_moist = 1005 + omega_in * 1860
    rho_air = psychrolib.GetMoistAirDensity(T_amb, omega_in, P_amb)

    print("density of air", rho_air)

    # Air properties
    mu_air = 1.8e-5
    air_speed = 2.0
    # m_da = 0.1 # fixed air mass flow
    m_da = np.round((rho_air * air_speed * 0.1),2) # mass flow of dry air based on velocity and area


    # Refrigerant side states
    T_evap_K = T_evap_air + 273.15
    T_cond_K = T_cond + 273.15

    P_evap = PropsSI("P", "T", T_evap_K, "Q", 1, refrigerant) # dryness fraction is used 1.0 for evaporater
    P_cond = PropsSI("P", "T", T_cond_K, "Q", 0, refrigerant)# dryness fraction is used 0.0 for condenser

    # Refrigerant states that do not depend on mass flow
    # State 1: saturated vapor at evaporator (outlet of evaporator)
    h1 = PropsSI("H", "T", T_evap_K, "Q", 1, refrigerant)
    s1 = PropsSI("S", "T", T_evap_K, "Q", 1, refrigerant)

    # State 3: saturated liquid at condenser outlet
    h3 = PropsSI("H", "T", T_cond_K, "Q", 0, refrigerant)
    
     # For throttling h4 = h3 (isenthalpic)
    h4 = h3

    PR = P_cond / P_evap
    eta = isentropic_eff(PR)
    h2s = PropsSI("H", "P", P_cond, "S", s1, refrigerant)
    h2 = h1 + (h2s - h1) / eta

    q_evap = h1 - h4
    mu_ref = PropsSI("V", "T", T_evap_K, "Q", 1, refrigerant)   # Dynamic viscosity of refrigerant at evaporator condition
    
    results = []



    # ---- Air-side required load ----

    # Reynolds number
    mu_air = PropsSI("V", "T", (T_amb+273), "P", P_amb, 'Air') 
    velocity = m_da / (rho_air * A_tube)
    Re = rho_air * velocity * D_tube / mu_air
    regime = "Laminar" if Re < 2300 else "Transitional" if Re < 4000 else "Turbulent"
    # --- Air side heat transfer (for future use) ---
    k_air = 0.026      # W/m-K
    Pr = 0.71

    if Re < 2300:
        Nu = 3.66
    else:
        Nu = 0.023 * (Re ** 0.8) * (Pr ** 0.4)

    h_air = Nu * k_air / D_tube


    cp_moist = 1005 + omega_in * 1860 # taken as dry air + water vapour 
    T_target = T_dew - 5

    Q_s = m_da * cp_moist * (T_amb - T_target)
    water = m_da * (omega_in - omega_out)
    Q_l = water * L_vap
    Q_required = Q_s + Q_l
    C_air = m_da * cp_moist
    NTU = U * A / C_air
    eps = 1 - np.exp(-NTU)

    results = []

    for m_ref in m_ref_vals:

        Q_ref = m_ref * q_evap
        Q_actual = eps * Q_ref

        W_comp = m_ref * (h2 - h1)
        COP = Q_actual / W_comp if W_comp > 0 else np.nan

        error = abs(Q_actual - Q_required)



        # Water production
        water_hr = water * 3600
        water_per_kWh = water_hr / (W_comp / 1000) if W_comp > 0 else np.nan

        # Reynolds number for refrigerant(inside evaporator tubes) ---
        Re_ref = m_ref * D_tube / (mu_ref * A_tube)
        regime_ref = "Laminar" if Re < 2300 else "Transitional" if Re < 4000 else "Turbulent"
        

        results.append([m_ref, m_da, Q_actual, Q_required, error, water_hr, water_per_kWh,
                            NTU, eps,W_comp,
                            COP, Re_ref, regime_ref])

        cols = ["m_ref", "m_da", "Q_W", "Q_required", "Error", "water_kg_hr", "water_kg_per_kWh",
                "NTU", "effectiveness", "W_comp",
                "COP", "Re", "flow_regime"]
    

    df = pd.DataFrame(results, columns=cols)
    out_name = "awg_results.csv"
    df.to_csv(out_name, index=False)
    print(f"Saved results to {out_name}")


    return df


# ---------------- MAIN EXECUTION ----------------
if __name__ == "__main__":
    df = run()

    # Find optimal refrigerant mass flow (minimum error)
    best = df.loc[df["Error"].idxmin()]

    print("\n===== OPTIMAL RESULT =====")
    print(f"Optimal refrigerant mass flow: {best['m_ref']:.4f} kg/s")
    print(f"Cooling delivered: {best['Q_W']:.2f} W")
    print(f"Compressor Power: {best['W_comp']:.2f} W")
    print(f"COP: {best['COP']:.2f}")
    print(f"Water production: {best['water_kg_hr']:.2f} kg/hr")
