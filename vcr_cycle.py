"""
vcr_cycle.py

What it does:
- Fetches ambient temperature, pressure, and relative humidity from OpenWeatherMap (or fallback to manual input).
- Computes psychrometric properties (dew point, humidity ratio, moist-air enthalpy).
- For a range of dry-air mass flow rates, computes:
    * air-side cooling Q_dot (including sensible + latent)
    * refrigerant cycle state points (R134a via CoolProp)
    * refrigerant mass flow required
    * compressor power and COP
    * water production rate (kg/hr) and water per kWh
- Writes a CSV with one row per mass flow and prints the recommended mass flow
Assumptions (editable in code):
- Evaporator refrigerant saturation temperature = (dew_point - 10°C)
- Condenser refrigerant saturation temperature = (T_ambient + 3°C)
- Evaporator treated as a heat exchanger with refrigerant phase-change side (C_r ~= 0)
- Isentropic compressor efficiency is estimated by a simple function (editable)
- Dry-air mass flow sweep default: 0.1 to 5.0 kg/s (per step)
"""

import os
import math
import numpy as np
import pandas as pd
import requests
from CoolProp.CoolProp import PropsSI
from dotenv import load_dotenv

import psychrolib
psychrolib.SetUnitSystem(psychrolib.SI)

# -------- USER / API CONFIG ----------
# Option A: use OpenWeatherMap (recommended)
# Get your API key from https://openweathermap.org/api (Current Weather Data)
load_dotenv()  # loads variables from .env
OWM_API_KEY = os.getenv("OWM_API_KEY") # not exposing api key directly as it's confidential.
if not OWM_API_KEY:
    print("⚠️ No OpenWeatherMap API key found in environment variables.")


# Provide location: either (lat, lon) tuple or city name string
# Examples:
# LOCATION = ("22.5726","88.3639")   # lat,lon for Kolkata
# LOCATION = "Kolkata,IN"
LOCATION = ("22.5726", "88.3639")

# If API not provided or fails, the script will ask for manual input (T_amb_C, P_amb_Pa, RH)
# ------------------------------------

# -------- PHYSICAL / DESIGN DEFAULTS (editable) ----------
refrigerant = "R134a"
# Temperatures offsets (°C) these are usually determined based on the model
EVAP_DEW_OFFSET = -10.0     # T_evap_air = dew_point + EVAP_DEW_OFFSET (note: offset is negative here)
COND_AMB_OFFSET = +3.0      # T_cond = T_ambient + COND_AMB_OFFSET

# Heat exchanger UA (if you want to compute NTU; not strictly used if we compute by enthalpy change)
UA_evap = 1200.0  # W/K (placeholder, not strongly used in current energy balance)
cp_dry_air = 1005.0 # J/kg-K 

# Isentropic efficiency estimate function (simple, editable)
def estimate_isentropic_efficiency(pressure_ratio):
    # baseline 0.75 for PR~2, reduce with PR increase (clipped)
    eta = 0.78 - 0.03 * max(0, (pressure_ratio - 2.0))
    return float(np.clip(eta, 0.60, 0.85))

# Sweep parameters (dry-air mass flow in kg/s) ranges of mass flow rates that are possibly optimal
m_dot_da_start = 0.1
m_dot_da_stop  = 5.0
m_dot_da_step  = 0.1

# Electricity cost (optional) for energy-per-hour metrics
ELECTRICITY_RATE_PER_KWH = 10.0  # INR per kWh (set if desired)

# ------------------------------------

# ---------- Helper: fetch weather ----------
def fetch_weather_openweathermap(location, api_key):
    """
    Accepts either "city,country" or (lat, lon) tuple of strings/floats.
    Returns a dict with 'T_ambient_C', 'P_amb_Pa', 'RH'
    """
    if api_key is None or api_key == "":
        raise ValueError("No OpenWeatherMap API key provided.")
    base = "https://api.openweathermap.org/data/2.5/weather"
    if isinstance(location, (list, tuple)):
        lat, lon = location
        params = {"lat": lat, "lon": lon, "appid": api_key, "units": "metric"}
    else:
        params = {"q": location, "appid": api_key, "units": "metric"}
    r = requests.get(base, params=params, timeout=10)
    r.raise_for_status()
    data = r.json()
    T_C = float(data["main"]["temp"])
    P_hpa = float(data["main"]["pressure"])  # hPa (or mbar)
    RH = float(data["main"]["humidity"]) / 100.0
    return {"T_ambient_C": T_C, "P_amb_Pa": P_hpa * 100.0, "RH": RH, "raw": data}

# ---------- Main ----------
def run_simulation(location=LOCATION, api_key=OWM_API_KEY):
    # 1) fetch ambient
    try:
        weather = fetch_weather_openweathermap(location, api_key)
        print("Fetched weather from OpenWeatherMap.")
    except Exception as e:
        print("OpenWeatherMap fetch failed:", e)
        # Fallback — ask user via input (non-interactive fallback values)
        # NOTE: If running non-interactively, edit these default fallbacks as needed.
        T_ambient_C = float(input("Enter ambient dry-bulb temperature (°C): "))
        P_amb_Pa = float(input("Enter ambient pressure (Pa): "))
        RH = float(input("Enter ambient relative humidity (0-1): "))
        weather = {"T_ambient_C": T_ambient_C, "P_amb_Pa": P_amb_Pa, "RH": RH}

    T_amb_C = weather["T_ambient_C"]
    P_amb_Pa = weather["P_amb_Pa"]
    RH = weather["RH"]

    # Derived psychrometric quantities
    Tdb = float(T_amb_C)
    T_dew = psychrolib.GetTDewPointFromRelHum(Tdb, RH)  # °C
    # evaporator air leaving temperature (assumed)
    T_evap_air_C = T_dew + EVAP_DEW_OFFSET
    T_cond_C = T_amb_C + COND_AMB_OFFSET

    print(f"Ambient: {T_amb_C:.2f} °C, RH: {RH*100:.1f} %, Dew point: {T_dew:.2f} °C")
    print(f"Assumed evaporator air temp (out): {T_evap_air_C:.2f} °C")
    print(f"Assumed condenser saturation temp: {T_cond_C:.2f} °C")

    # Precompute enthalpies for moist air inlet/outlet at a given humidity ratio
    # But humidity ratio depends on mass flow; we'll compute inside loop (as omega_in depends on Tdb,RH)
    # However omega_in is independent of m_dot
    omega_in = psychrolib.GetHumRatioFromRelHum(Tdb, RH, P_amb_Pa)  # kg water/kg dry air
    omega_sat_evap = psychrolib.GetHumRatioFromRelHum(T_evap_air_C, 1.0, P_amb_Pa)
    h_in = psychrolib.GetMoistAirEnthalpy(Tdb, omega_in)   # J/kg dry air
    h_out_saturated = psychrolib.GetMoistAirEnthalpy(T_evap_air_C, omega_sat_evap)

    # Refrigerant-side temps (Kelvin)
    T_evap_R_K = T_evap_air_C + 273.15
    T_cond_R_K = T_cond_C + 273.15

    # Refrigerant pressures (saturation)
    P_evap = PropsSI("P", "T", T_evap_R_K, "Q", 1.0, refrigerant) # dryness fraction is used 1.0 for evaporater
    P_cond = PropsSI("P", "T", T_cond_R_K, "Q", 0.0, refrigerant) # dryness fraction is used 0.0 for condenser

    # Refrigerant states that do not depend on mass flow
    # State 1: saturated vapor at evaporator (outlet of evaporator)
    h1 = PropsSI("H", "T", T_evap_R_K, "Q", 1.0, refrigerant)  # J/kg
    s1 = PropsSI("S", "T", T_evap_R_K, "Q", 1.0, refrigerant)

    # State 3: saturated liquid at condenser outlet
    h3 = PropsSI("H", "T", T_cond_R_K, "Q", 0.0, refrigerant)

    # For throttling h4 = h3 (isenthalpic)
    h4 = h3

    # Preparation of sweep
    # Preparation of sweep (rounded to avoid floating-point drift)
    m_vals = np.round( np.arange(m_dot_da_start, m_dot_da_stop, m_dot_da_step), 2)
    results = []

    for m_da in m_vals:
        # Air-side cooling load (per second) using moist-air enthalpy (J/s)
        Q_dot_air_W = m_da * (h_in - h_out_saturated)  # W (J/s)
        if Q_dot_air_W <= 0:
            # no cooling possible (e.g., input selection), skip
            results.append({
                "m_dot_da_kg_s": m_da,
                "Q_dot_air_W": Q_dot_air_W,
                "m_ref_kg_s": 0.0,
                "W_comp_W": 0.0,
                "COP": 0.0,
                "water_kg_hr": max(0.0, m_da * (omega_in - omega_sat_evap) * 3600.0),
                "water_kg_per_kWh": 0.0
            })
            continue

        # Refrigerant: compute isentropic compression to condenser pressure
        # h2s at (P_cond, s1)
        try:
            h2s = PropsSI("H", "P", P_cond, "S", s1, refrigerant)
        except Exception as e:
            # if property call fails, set NaNs and continue
            h2s = float("nan")

        # pressure ratio
        PR = float(P_cond / P_evap) if P_evap > 0 else 1.0
        eta_isentropic = estimate_isentropic_efficiency(PR)
        # actual h2
        h2 = h1 + (h2s - h1) / eta_isentropic

        # Specific evaporator refrigeration per kg refrigerant (J/kg)
        q_evap_specific = h1 - h4
        if q_evap_specific <= 0:
            m_ref = float("nan")
            W_comp = float("nan")
            COP = float("nan")
        else:
            # refrigerant mass flow required to absorb Q_dot_air
            m_ref = Q_dot_air_W / q_evap_specific  # kg refrigerant / s

            # compressor power total (W): m_ref * (h2 - h1)
            W_comp = m_ref * (h2 - h1)

            # COP
            COP = Q_dot_air_W / W_comp if W_comp > 0 else float("nan")

        # water collection rate (kg/hr) (per kg dry air flow)
        water_kg_hr = max(0.0, m_da * (omega_in - omega_sat_evap) * 3600.0)

        # water per kWh (kg produced divided by kWh consumed)
        energy_kW = W_comp / 1000.0 if W_comp > 0 else float("nan")
        water_kg_per_kWh = (water_kg_hr / energy_kW) if (energy_kW and energy_kW > 0) else float("nan")

        results.append({
            "m_dot_da_kg_s": m_da,
            "Q_dot_air_W": Q_dot_air_W,
            "omega_in": omega_in,
            "omega_out_sat": omega_sat_evap,
            "h_air_in_J_per_kg_dry": h_in,
            "h_air_out_J_per_kg_dry": h_out_saturated,
            "T_evap_air_C": T_evap_air_C,
            "T_cond_C": T_cond_C,
            "P_evap_Pa": P_evap,
            "P_cond_Pa": P_cond,
            "h1_Jkg": h1,
            "h2s_Jkg": h2s,
            "h2_Jkg": h2,
            "h3_Jkg": h3,
            "h4_Jkg": h4,
            "q_evap_specific_J_per_kg_ref": q_evap_specific,
            "m_ref_kg_s": m_ref,
            "W_comp_W": W_comp,
            "COP": COP,
            "water_kg_hr": water_kg_hr,
            "water_kg_per_kWh": water_kg_per_kWh,
            "eta_isentropic": eta_isentropic
        })

    df = pd.DataFrame(results)
    # Save CSV
    out_name = "air_water_refrigerator_sweep_results.csv"
    df.to_csv(out_name, index=False, float_format="%.4f")
    print(f"Saved results to {out_name}")

    # Pick best mass flow by water per kWh (descending)
    # ignore NaNs
    df_valid = df[df["water_kg_per_kWh"].notnull() & (df["water_kg_per_kWh"] > 0)]
    if not df_valid.empty:
        best_row = df_valid.sort_values("water_kg_per_kWh", ascending=False).iloc[0]
        print("Best mass flow by water/kg per kWh:")
        print(best_row[["m_dot_da_kg_s", "water_kg_hr", "water_kg_per_kWh", "COP", "W_comp_W", "Q_dot_air_W"]])
    else:
        print("No valid results to determine best mass flow.")

    return df

if __name__ == "__main__":
    df_results = run_simulation()
    # quick summary head
    print(df_results.head(8))
    print("Done.")
