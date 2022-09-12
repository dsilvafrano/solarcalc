import numpy as np
import pandas as pd
from scipy.optimize import minimize, basinhopping
import sqlite3
import logger
import numpy_financial as npf
import time
# from shgo import shgo
import requests
from datetime import date, datetime, timedelta
from calendar import monthrange

# from WhatIfAnalysis import GoalSeek
import json

#Inputs from users
sload = 20 # from user
tariff_id = 1 # consumer type
voltage_id = 1 # voltage type
voltage = "LT"
tariff = 'Domestic'
residence_type = 'Independent House'
metering_id = 7
state_id = 1
state = "Tamil Nadu"
metering_type = "Net Metering"
weekend_consumption_change = 0
weekend_consumption_separate = 0 # if this is 1, it means there is weekend consumption
load_input_type = "average_monthly"
avg_monthly = 5000
# in case of monthwise load_input_type
mc1 = 5000
mc2 = 5000
mc3 = 5000
mc4 = 5000
mc5 = 5000
mc6 = 5000
mc7 = 5000
mc8 = 5000
mc9 = 5000
mc10 = 5000
mc11 = 5000
mc12 = 5000

tou_select = 2
nyr = 26
solar = True
battery = False
# x1=np.zeros(2, dtype=float)
# x1[0] = 10 # user input solar capacity
# x1[1] = 10 # user input storage capacity
der_deg = 0.01  # solar degradation
bat_type = 1  # battery type = 1 fo Li ion, 0 for lead acid
socin = 0.6
socmax = 0.9
#Ma
batstatus = pd.read_csv("dispatch2.csv", header=None)  # need to be replaced with the actual dispatch strategy
# print(batstatus[24:48])

solarpv_subsidy = 0

# starting time
start = time.time()

# fetch solar data using latitude and longitude (get the latitude and longitude using the pincode data
latitude = 12.0055
longitude = 79.8089
response = requests.get("https://developer.nrel.gov/api/pvwatts/v6.json?api_key"
                        "=t8EOwc4wOkyTPfitOHDGgaPqUsWYVw9B3JxaG6Fu&lat=" + str(latitude) + "&lon=" + str(
    longitude) + "&system_capacity=1"
                 "&azimuth=180&tilt=12&array_type=1&module_type=1&losses=11&timeframe=hourly")
outpts = response.json()['outputs']
acpwr = outpts['ac']
# print(acpwr)
solarp = pd.DataFrame()
solarp['load'] = acpwr
# print(solarp[0:10])
# print(solarp.columns)
solarp = solarp / 1000
# print(solarp[0:10])



# create connection to DB
random_no = 0.3555
def create_connection():
    conn = None
    try:
        conn = sqlite3.connect("db6_updated.sqlite3")
    except ValueError as e:
        logger.debug(random_no + " " + e)
    return conn

conn = create_connection()

# fetch the dc ac ratio from database
dc_ac_ratio_q = pd.read_sql_query(
    "select value from assumptions_pvwatts where parameter = 'dc_ac_ratio'", conn)
pysam_dc_ac_ratio = float(dc_ac_ratio_q.values[0])
# fetch the inverter replacement year from databas
inv_replace_year_q = pd.read_sql_query(
    "select value from assumptions_pvwatts where parameter = 'inverter_replacement_year'", conn)
inv_replace_year = int(inv_replace_year_q.values[0])

# find the battery replacement year based on the lifecycle of batteries
if bat_type == 1:
    battery_cycle_q = pd.read_sql_query(
        "select value from assumptions_battwatts where parameter = 'battery_cyclelife_liion'", conn)

    battery_cycle_v = int(battery_cycle_q.values[0])
    battery_dod_q = pd.read_sql_query(
        "select value from assumptions_battwatts where parameter = 'battery_dod_liion'", conn)
    dod = int(battery_dod_q.values[0])/100
elif bat_type == 0:
    battery_cycle_q = pd.read_sql_query(
        "select value from assumptions_battwatts where parameter = 'battery_cyclelife_pbacid'", conn)
    battery_cycle_v = int(battery_cycle_q.values[0])
    battery_dod_q = pd.read_sql_query(
        "select value from assumptions_battwatts where parameter = 'battery_dod_pbacid'", conn)
    dod = int(battery_dod_q.values[0])/100
print('DOD :', dod)

socmin = 1 - dod
print('SOC min', socmin)
battery_replace_year = int(battery_cycle_v / 365)
numb_battery_replacement = int(nyr / (battery_replace_year + 1))
rep_yrs_battery = [i * (battery_replace_year + 1) for i in range(1, numb_battery_replacement + 1)]
rep_yrs_inverter = [inv_replace_year]

# loan related data

# debt fracton fetching from database
debt_fraction_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'debt_fraction'", conn)
# logger.debug("pysam_debt_fraction ="+str((debt_fraction_q.values[0])))
pysam_debt_fraction = float(debt_fraction_q.values[0])

# Loan Rate from database
loan_rate_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'loan_rate'", conn)
# logger.debug(random_no+" "+"loan_rate_q ="+str((loan_rate_q.values[0])))
loan_rate = float(loan_rate_q.values[0]) / 100

# loan period
loan_period_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'loan_term'", conn)
# logger.debug(random_no+" "+"loan_rate_q ="+str((loan_rate_q.values[0])))
loan_period = float(loan_period_q.values[0])

# cost escalation/inflation rate
cost_esc_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'inflation_rate'", conn)
# logger.debug(random_no+" "+"cost_esc_q ="+str((lcost_esc_q.values[0])))
cost_esc = cost_esc_q.values[0] / 100
print('Cost escalation :', cost_esc)

# get discount rate for calculating npv
dis_factor_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'real_discount_rate_show'",
                                 conn)
# logger.debug(random_no+" "+"cost_esc_q ="+str((lcost_esc_q.values[0])))
dis_factor = dis_factor_q.values[0] / 100
print('Discount rate :', dis_factor)

# cost escalation/inflation rate
cost_esc_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'inflation_rate'", conn)
# logger.debug(random_no+" "+"cost_esc_q ="+str((lcost_esc_q.values[0])))
cost_esc = cost_esc_q.values[0] / 100
print('Cost escalation :', cost_esc)

# load escalation
load_esc_q = pd.read_sql_query("select value from assumptions_grid where parameter = 'load_escalation'", conn)
# logger.debug(random_no+" "+"cost_esc_q ="+str((lcost_esc_q.values[0])))
load_esc = load_esc_q.values[0] / 100
print('Load escalation:', load_esc)

# get electricity charges and fixed charges
# cost_energy_q = pd.read_sql_query(
#     "select energy_charge from optimizer_charges where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(
#         voltage_id) + " and state_id=" + str(state_id), conn)
# cost_energy = float(cost_energy_q.values[0])
# print("cost_energy =", str(cost_energy))

# Build a matrix for energy cost for tiered calculation
slab_id_q = pd.read_sql_query("select slab_id, period, tier, min, maximum, energy_charge from slabs_mapping where state_id =" + str(state_id) + " and tarriff_type_id = " + str(tariff_id) + " and voltage_type_id = " + str(voltage_id) + " and metering_type_id = " + str(metering_id),conn)
slab_id_t = slab_id_q
# print('The EC table is:', slab_id_t)
slab_id_m = slab_id_t['slab_id'].max()
# print('The max slab is', slab_id_m)
tier_m = slab_id_t['tier'].max()
# print('The tier slab is', tier_m)


# get fixed charge calculation type
charge_calculation_q = pd.read_sql_query("select charge_calculation from fixedcharge where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
charge_calculation = int(charge_calculation_q.values[0])
if int(charge_calculation) == 0:
    fixed_charge = 0
elif charge_calculation == 1:
    fixed_charges_q = pd.read_sql_query("select fixed_charge from fixedcharge where state_id =" + str(state_id) + " and voltage_id = " + str(voltage_id),conn)
    fixed_charge = float(fixed_charges_q.values[0])
elif charge_calculation == 2:

#Tier check for tiers larger than 1
    tier_check_q = pd.read_sql_query("select tier from fixedcharge where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
    tier_check = tier_check_q
    # print('Tier check:', tier_check)
    max_tier_check = tier_check['tier'].max()
    # print('Max tier check:', max_tier_check)
    min_tier_check = tier_check['tier'].min()
    # print('Min tier check:', min_tier_check)
# Fetch fixed charges from database & Tiered calculation
    if max_tier_check == 1:
        fixed_charges_q = pd.read_sql_query("select fixed_charge from fixedcharge where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
        fixed_charge = float(fixed_charges_q.values[0])
        fixed_charge = fixed_charge * sload
    else:
        sl_check_q = pd.read_sql_query("select fixed_charge, sanctioned_load from fixedcharge where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
        sl_fc_check = sl_check_q
        # print('SL:', sl_fc_check)
        # print('FC:', sl_fc_check['fixed_charge'][1])
        # print('SL:', sl_fc_check['sanctioned_load'][1])

        fixed_charge = 0
        fc_sl = sload
        fc_temp = 0
        for i in range (0, len(sl_fc_check)):
            # print('new sl value:', fc_sl)
            if fc_sl >= float(sl_fc_check['sanctioned_load'][i]):
                fc_temp = fc_temp + (sl_fc_check['fixed_charge'][i] * sl_fc_check['sanctioned_load'][i])
                fc_sl = fc_sl - float(sl_fc_check['sanctioned_load'][i])
                fc_sl = fc_sl
                # print('FC temp:', fc_sl)
            else:
                fc_temp = fc_temp + (sl_fc_check['fixed_charge'][i] * fc_sl)
                fc_sl = fc_sl - fc_sl
            # print('FC value:', fc_temp)
        fixed_charge = fc_temp
        fixed_charge = fixed_charge

fixed_charge_a = fixed_charge*12
fixed_charge_h = fixed_charge_a/8760

# print("fixed charge =", str(fixed_charge))
# print("fixed charge annual =", str(fixed_charge_a))
# print("fixed charge per hour =", str(fixed_charge_h))

# Selection on network charge and compensation rate from network charge table in SQL
def network_charge_fetch(user_pv_capacity):
    if state_id == 1 and voltage_id == 1:
        if float(user_pv_capacity) <= 10:
            network_charge_df = pd.read_sql_query("select network_charge from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(10), conn)
            comp_charge_df = pd.read_sql_query("select compensation_rate from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(10), conn)

        elif 10 < float(user_pv_capacity) <= 112:
            network_charge_df = pd.read_sql_query("select network_charge from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(112), conn)
            comp_charge_df = pd.read_sql_query("select compensation_rate from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(112), conn)

        else:
            network_charge_df = pd.read_sql_query("select network_charge from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(999), conn)

            comp_charge_df = pd.read_sql_query("select compensation_rate from network_charge where state_id=" + str(
                state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
                tariff_id) + " and voltage_id=" + str(voltage_id) + " and max_pv=" + str(999), conn)
    else:
        network_charge_df = pd.read_sql_query("select network_charge from network_charge where state_id=" + str(
            state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
            tariff_id) + " and voltage_id=" + str(voltage_id), conn)
        comp_charge_df = pd.read_sql_query("select compensation_rate from network_charge where state_id=" + str(
            state_id) + " and metering_type_id=" + str(metering_id) + " and tariff_id=" + str(
            tariff_id) + " and voltage_id=" + str(voltage_id), conn)

    network_charges = float(network_charge_df.values[0])
    compensation_rates = float(comp_charge_df.values[0])
    # print('compensation rate:', compensation_rates)
    # print('network_charges:', network_charges)

    return network_charges, compensation_rates

# Fetching financial for solar, inverter & battery from cost suggestion table in SQL
def financial_fetch(user_pv_capacity):

    df_for_state_id = pd.read_sql_query("select id from state WHERE name=" + f'"{state}"', conn)
    state_id = df_for_state_id['id'].values[0]

    df_for_consumer_id = pd.read_sql_query("select id from tariff_type WHERE name=" + f'"{tariff}"', conn)
    consumer_id = df_for_consumer_id['id'].values[0]

    df_for_voltage_id = pd.read_sql_query("select id from voltage_type WHERE name=" + f'"{voltage}"', conn)
    voltage_id = df_for_voltage_id['id'].values[0]

    df_for_metering_id = pd.read_sql_query(
        "select id from metering_type WHERE name=" + f'"{metering_type}"' + "and state_id = " + f'"{str(state_id)}"', conn);
    metering_id = df_for_metering_id['id'].values[0]

    # if df_for_max_pv.shape[0] == 0:
    df_max_pv = pd.read_sql_query(
        "select max_pv,min_pv,network_charge,compensation_rate from network_charge WHERE state_id=" + f'"{state_id}"' + " and tariff_id=" + f'"{consumer_id}"' + " and voltage_id=" + f'"{voltage_id}"' + " and metering_type_id=" + f'"{metering_id}"',
        conn)
    max_limit_pv = df_max_pv['max_pv'].max()
    min_limit_pv = df_max_pv['min_pv'].min()

    if max_limit_pv > float(sload):
        max_limit_pv = sload
    if str(metering_id) == "8":
        if min_limit_pv < float(sload): #This is as SL for TN is 112kW and min capacity under Gross metering is 150kW
            min_limit_pv = 150


    df_for_cost_suggestion = pd.read_sql_query(
        "select capital_cost,inverter,pv_cost,hybrid_inverter,li_ion,lead_acid from cost_suggestion WHERE " + str(sload) + " between min_pv_cap and max_pv_cap",
        conn)


    max_limit_pv = str(max_limit_pv)
    min_limit_pv = str(min_limit_pv)
    capital_cost = str(df_for_cost_suggestion['capital_cost'].values[0])
    inverter_cost = str(df_for_cost_suggestion['inverter'].values[0])
    hybrid_inverter = str(df_for_cost_suggestion['hybrid_inverter'].values[0])
    pv_cost = str(df_for_cost_suggestion['pv_cost'].values[0])
    li_ion = str(df_for_cost_suggestion['li_ion'].values[0])
    lead_acid = str(df_for_cost_suggestion['lead_acid'].values[0])

    print('The capital cost is:', capital_cost)

    return max_limit_pv, min_limit_pv, capital_cost, inverter_cost, hybrid_inverter, pv_cost, li_ion, lead_acid
    #print(
        # "select max_pv from network_charge WHERE state_id=" + f'"{state_id}"' + " and tariff_id=" + f'"{consumer_id}"' + " and voltage_id=" + f'"{voltage_id}"' + " and metering_type_id=" + f'"{metering_id}"')


    #max_limit_pv = d.read_sql_query("select max_pv from network_charge WHERE state_id=" + f'"{state_id}"' + " and tariff_id=" + f'"{consumer_id}"' + " and voltage_id=" + f'"{voltage_id}"' + " and metering_type_id=" + f'"{metering_id}"')

    #min_limit_pv = d.read_sql_query("select max_pv from network_charge WHERE state_id=" + f'"{state_id}"' + " and tariff_id=" + f'"{consumer_id}"' + " and voltage_id=" + f'"{voltage_id}"' + " and metering_type_id=" + f'"{metering_id}"')

# Calculation of initial investment
def investmentcost_calculate(system_capacity, bat_sim_kwh):
    print('solar size:', system_capacity)
    if solar & battery:
        if bat_type == 1:
            if tariff == 'Domestic' and residence_type == 'Independent House' and system_capacity < 500:
                total_installation_cost = (float(system_capacity) * float(financial_fetch(sload)[5]) + float(bat_sim_kwh) * int(
                    financial_fetch(sload)[6]) + int(financial_fetch(sload)[4]) * float(system_capacity) / float(pysam_dc_ac_ratio)) - int(
                    solarpv_subsidy) * system_capacity
            else:
                total_installation_cost = (
                        (float(system_capacity) * float(financial_fetch(sload)[5])) + (float(bat_sim_kwh) * int(financial_fetch(sload)[6])) +
                        (int(financial_fetch(sload)[4]) * (float(system_capacity) / float(pysam_dc_ac_ratio))))
        else:
            if tariff == 'Domestic' and residence_type == 'Independent House' and system_capacity < 500:
                total_installation_cost = (float(system_capacity) * float(financial_fetch(sload)[5]) + float(bat_sim_kwh) * int(
                    financial_fetch(sload)[7]) + int(financial_fetch(sload)[4]) * float(system_capacity) / float(pysam_dc_ac_ratio)) - int(
                    solarpv_subsidy) * system_capacity
            else:
                total_installation_cost = (
                        (float(system_capacity) * float(financial_fetch(sload)[5])) + (float(bat_sim_kwh) * int(financial_fetch(sload)[7])) +
                        (int(financial_fetch(sload)[4]) * (float(system_capacity) / float(pysam_dc_ac_ratio))))
    else:
        if tariff == 'Domestic' and residence_type == 'Independent House' and system_capacity < 500:
            total_installation_cost = (float(system_capacity) * float(financial_fetch(sload)[2])) - int(
                solarpv_subsidy) * system_capacity
        else:
            total_installation_cost = (float(system_capacity) * float(financial_fetch(sload)[2]))

    return total_installation_cost


# find the replacement cost of inverter and battery
def replacement_cost(system_capacity, bat_sim_kwh):
    rep_battery_cost = np.zeros(nyr)
    rep_inverter_cost = np.zeros(nyr)
    if bat_type == 1:
        battery_cost = float(financial_fetch(sload)[6])
    else:
        battery_cost = float(financial_fetch(sload)[7])

    if solar & battery:
        for k in range(0, nyr):
            i = k + 1
            if i in rep_yrs_battery:
                rep_battery_cost[k] = battery_cost * bat_sim_kwh
            if i in rep_yrs_inverter:
                rep_inverter_cost[k] = float(financial_fetch(sload)[3]) * system_capacity
    else:
        for k in range(0, nyr):
            i = k + 1
            if i in rep_yrs_inverter:
                rep_inverter_cost[k] = float(financial_fetch(sload)[3]) * system_capacity
    return rep_battery_cost, rep_inverter_cost




# Inputs from Solsavi
user_load = np.array(([0] * 8760), dtype=float)

n = len(user_load)  # number of load values

# Define load distribution for 24hrs weekday
weekday_consumption_6to10 = 20
weekday_consumption_10to18 = 40
weekday_consumption_18to22 = 30
weekday_consumption_22to6 = 10


# Define load distribution for 24hrs weekend


weekend_consumption_6to10 = 20
weekend_consumption_10to18 = 20
weekend_consumption_18to22 = 40
weekend_consumption_22to6 = 20

weekday_consumption_6to10n = weekday_consumption_6to10 / ((
                weekday_consumption_6to10 + weekday_consumption_10to18 + weekday_consumption_18to22 + weekday_consumption_22to6))
weekday_consumption_10to18n = weekday_consumption_10to18 / ( (
                weekday_consumption_6to10 + weekday_consumption_10to18 + weekday_consumption_18to22 + weekday_consumption_22to6))
weekday_consumption_18to22n = weekday_consumption_18to22 / ( (
                weekday_consumption_6to10 + weekday_consumption_10to18 + weekday_consumption_18to22 + weekday_consumption_22to6))
weekday_consumption_22to6n = weekday_consumption_22to6 / ( (
                weekday_consumption_6to10 + weekday_consumption_10to18 + weekday_consumption_18to22 + weekday_consumption_22to6))

if int(weekend_consumption_separate) == 1:
    weekend_consumption_6to10n = weekend_consumption_6to10 / ((
            weekend_consumption_6to10 + weekend_consumption_10to18 + weekend_consumption_18to22 + weekend_consumption_22to6))
    weekend_consumption_10to18n = weekend_consumption_10to18 / ((
            weekend_consumption_6to10 + weekend_consumption_10to18 + weekend_consumption_18to22 + weekend_consumption_22to6))
    weekend_consumption_18to22n = weekend_consumption_18to22 / ((
            weekend_consumption_6to10 + weekend_consumption_10to18 + weekend_consumption_18to22 + weekend_consumption_22to6))
    weekend_consumption_22to6n = weekend_consumption_22to6 / ((
            weekend_consumption_6to10 + weekend_consumption_10to18 + weekend_consumption_18to22 + weekend_consumption_22to6))

#Populating the days of the month
def last_day_of_month(date_value):
    return date_value.replace(day = monthrange(date_value.year, date_value.month)[1])

# In the case where load_input_type = "Average monthly"
# Populate data in weekday
if load_input_type == "average_monthly":
    if int(weekend_consumption_separate) == 0:
        user_load[6:10] = round((weekday_consumption_6to10n * avg_monthly * 12 / (365 * 4)), 3)
        user_load[10:18] = round((weekday_consumption_10to18n * avg_monthly * 12 / (365 * 8)), 3)
        user_load[18:22] = round((weekday_consumption_18to22n * avg_monthly * 12 / (365 * 4)), 3)
        user_load[22:24] = round((weekday_consumption_22to6n * avg_monthly * 12 * 0.25 / (365 * 2)), 3)
        user_load[0:6] = round((weekday_consumption_22to6n * avg_monthly * 12 * 0.75 / (365 * 6)), 3)
        a = 0
        b = 24
        for i in range(1, 365):
            x = a
            y = b
            a = b
            b = b + 24
            user_load[a:b] = user_load[x:y]

        print("user load full matrix ",user_load)
        print(len(user_load))
    else: # Populate data in weekend
        weekday_24 = np.array(([0] * 24), dtype=float)
        weekend_24 = np.array(([0] * 24), dtype=float)

        weekday_daily = (avg_monthly * 12) / (261 + (104 * (weekend_consumption_change + 1)))
        weekend_daily = (weekend_consumption_change + 1) * weekday_daily

        # print(weekday_daily)
        # print(weekend_daily)

        # print(weekday_consumption_6to10n)

        weekday_24[6:10] = round((weekday_consumption_6to10n * (weekday_daily) / 4), 3)
        weekday_24[10:18] = round((weekday_consumption_10to18n * (weekday_daily) / 8), 3)
        weekday_24[18:22] = round((weekday_consumption_18to22n * (weekday_daily) / 4), 3)
        weekday_24[22:24] = round((weekday_consumption_22to6n * (weekday_daily * 0.25) / 2), 3)
        weekday_24[0:6] = round((weekday_consumption_22to6n * (weekday_daily * 0.75) / 6), 3)

        weekend_24[6:10] = round((weekend_consumption_6to10n * (weekend_daily) / 4), 3)
        weekend_24[10:18] = round((weekend_consumption_10to18n * (weekend_daily) / 8), 3)
        weekend_24[18:22] = round((weekend_consumption_18to22n * (weekend_daily) / 4), 3)
        weekend_24[22:24] = round((weekend_consumption_22to6n * (weekend_daily * 0.25) / 2), 3)
        weekend_24[0:6] = round((weekend_consumption_22to6n * (weekend_daily * 0.75) / 6), 3)

        # print(weekday_24[0:24])
        # print(weekend_24[0:24])

        year = date.today().year

        begin_year = date(year, 1, 1)
        end_year = date(year, 12, 31)
        one_day = timedelta(days=1)

        print(begin_year, end_year, one_day)

        a = 0
        b = 24

        next_day = begin_year
        for day in range(0, 365):  # includes potential leap year
            if next_day > end_year:
                break
            if (next_day.weekday() > 4):
                user_load[a:b] = weekend_24[0:24]
                a = b
                b = b + 24
            else:
                user_load[a:b] = weekday_24[0:24]
                a = b
                b = b + 24
            next_day += one_day
        # print("user load full matrix ", user_load)
        # print(len(user_load))
    # Populate yearly load to monthly
    Year1_monthly = []
    jan_1 = 0
    feb_1 = 0
    mar_1 = 0
    apr_1 = 0
    may_1 = 0
    jun_1 = 0
    jul_1 = 0
    aug_1 = 0
    sep_1 = 0
    oct_1 = 0
    nov_1 = 0
    dec_1 = 0

    for i in range(0, 8760):
        if i in range(0, 744):
            jan_1 = jan_1 + user_load[i]
        elif i in range(745,1416):
            feb_1 = feb_1 + user_load[i]
        elif i in range(1417,2160):
            mar_1 = mar_1 + user_load[i]
        elif i in range(2161,2880):
            apr_1 = apr_1 + user_load[i]
        elif i in range(2881,3624):
            may_1 = may_1 + user_load[i]
        elif i in range(3625,4344):
            jun_1 = jun_1 + user_load[i]
        elif i in range(4345,5088):
            jul_1 = jul_1 + user_load[i]
        elif i in range(5089,5832):
            aug_1 = aug_1 + user_load[i]
        elif i in range(5832,6552):
            sep_1 = sep_1 + user_load[i]
        elif i in range(6553,7296):
            oct_1 = oct_1 + user_load[i]
        elif i in range(7297,8016):
            nov_1 = nov_1 + user_load[i]
        else:
            dec_1  = dec_1 + user_load[i]

    Year1_monthly.append(jan_1)
    Year1_monthly.append(feb_1)
    Year1_monthly.append(mar_1)
    Year1_monthly.append(apr_1)
    Year1_monthly.append(may_1)
    Year1_monthly.append(jun_1)
    Year1_monthly.append(jul_1)
    Year1_monthly.append(aug_1)
    Year1_monthly.append(sep_1)
    Year1_monthly.append(oct_1)
    Year1_monthly.append(nov_1)
    Year1_monthly.append(dec_1)

    print(Year1_monthly)
    print(len(Year1_monthly))

elif load_input_type == "month_wise":
    mc = []
    mc.append(mc1)
    mc.append(mc2)
    mc.append(mc3)
    mc.append(mc4)
    mc.append(mc5)
    mc.append(mc6)
    mc.append(mc7)
    mc.append(mc8)
    mc.append(mc9)
    mc.append(mc10)
    mc.append(mc11)
    mc.append(mc12)
    year = date.today().year;
    if int(weekend_consumption_separate) == 0:
        a = 6
        b = 10
        c = 18
        d = 22
        e = 24
        f = 0
        for i in range(1, 13):
            # print("Completed Month ", i)
            days_in_month = monthrange(2018, i)[1]
            for j in range(1, days_in_month + 1):
                user_load[a:b] = round((weekday_consumption_6to10n * mc[i - 1] / (days_in_month * 4)), 3)
                user_load[b:c] = round((weekday_consumption_10to18n * mc[i - 1] / (days_in_month * 8)), 3)
                user_load[c:d] = round((weekday_consumption_18to22n * mc[i - 1] / (days_in_month * 4)), 3)
                user_load[d:e] = round((weekday_consumption_22to6n * mc[i - 1] * 0.25 / (days_in_month * 2)), 3)
                user_load[f:a] = round((weekday_consumption_22to6n * mc[i - 1] * 0.75 / (days_in_month * 6)), 3)
                a = a + 24
                b = b + 24
                c = c + 24
                d = d + 24
                e = e + 24
                f = f + 24
    else:
        a = 6
        b = 10
        c = 18
        d = 22
        e = 24
        f = 0
        for i in range(1, 13):
            # print("Completed Month ", i)

            days_in_month = monthrange(year, i)[1]
            given_date = datetime(year=year, month=i, day=1).date()
            first_day_of_month = given_date.replace(day=1)
            weekdays_in_month = np.busday_count(first_day_of_month, last_day_of_month(given_date))
            weekends_in_month = days_in_month - weekdays_in_month

            # print("weekdays ", weekdays_in_month)
            # print("weekends ", weekends_in_month)

            consumption_weekday = mc[i - 1] / (
                        weekdays_in_month + (weekends_in_month * (weekend_consumption_change + 1)))
            consumption_weekend = ((weekend_consumption_change + 1) * consumption_weekday)
            # print('consumption_weekday ', consumption_weekday)

            # print('consumption_weekend ', consumption_weekend)

            # print(weekend_consumption_6to10n)
            # print(weekend_consumption_10to18n)
            # print(weekend_consumption_18to22n)
            # print(weekend_consumption_22to6n)

            weekday_24 = np.array(([0] * 24), dtype=float)
            weekend_24 = np.array(([0] * 24), dtype=float)

            weekday_24[6:10] = round((weekday_consumption_6to10n * consumption_weekday / 4), 3)
            weekday_24[10:18] = round((weekday_consumption_10to18n * consumption_weekday / 8), 3)
            weekday_24[18:22] = round((weekday_consumption_18to22n * consumption_weekday / 4), 3)
            weekday_24[22:24] = round((weekday_consumption_22to6n * consumption_weekday * 0.25 / 2), 3)
            weekday_24[0:6] = round((weekday_consumption_22to6n * consumption_weekday * 0.75 / 6), 3)

            weekend_24[6:10] = round((weekend_consumption_6to10n * consumption_weekend / 4), 3)
            weekend_24[10:18] = round((weekend_consumption_10to18n * consumption_weekend / 8), 3)
            weekend_24[18:22] = round((weekend_consumption_18to22n * consumption_weekend / 4), 3)
            weekend_24[22:24] = round((weekend_consumption_22to6n * consumption_weekend * 0.25 / 2), 3)
            weekend_24[0:6] = round((weekend_consumption_22to6n * consumption_weekend * 0.75 / 6), 3)

            year = date.today().year

            begin_year = date(year, 1, 1)
            end_year = date(year, 12, 31)
            one_day = timedelta(days=1)

            # load_value = [];

            # print(begin_year, end_year, one_day)

            a = 0
            b = 24

            next_day = begin_year
            for day in range(0, 365):  # includes potential leap year
                if next_day > end_year:
                    break
                if (next_day.weekday() > 4):
                    user_load[a:b] = weekend_24[0:24]
                    a = b
                    b = b + 24
                else:
                    user_load[a:b] = weekday_24[0:24]
                    a = b
                    b = b + 24
                next_day += one_day
    # Populate yearly load to monthly
    Year1_monthly = []
    jan_1 = 0
    feb_1 = 0
    mar_1 = 0
    apr_1 = 0
    may_1 = 0
    jun_1 = 0
    jul_1 = 0
    aug_1 = 0
    sep_1 = 0
    oct_1 = 0
    nov_1 = 0
    dec_1 = 0

    for i in range(0, 8760):
        if i in range(0, 744):
            jan_1 = jan_1 + user_load[i]
        elif i in range(745, 1416):
            feb_1 = feb_1 + user_load[i]
        elif i in range(1417, 2160):
            mar_1 = mar_1 + user_load[i]
        elif i in range(2161, 2880):
            apr_1 = apr_1 + user_load[i]
        elif i in range(2881, 3624):
            may_1 = may_1 + user_load[i]
        elif i in range(3625, 4344):
            jun_1 = jun_1 + user_load[i]
        elif i in range(4345, 5088):
            jul_1 = jul_1 + user_load[i]
        elif i in range(5089, 5832):
            aug_1 = aug_1 + user_load[i]
        elif i in range(5832, 6552):
            sep_1 = sep_1 + user_load[i]
        elif i in range(6553, 7296):
            oct_1 = oct_1 + user_load[i]
        elif i in range(7297, 8016):
            nov_1 = nov_1 + user_load[i]
        else:
            dec_1 = dec_1 + user_load[i]

    Year1_monthly.append(jan_1)
    Year1_monthly.append(feb_1)
    Year1_monthly.append(mar_1)
    Year1_monthly.append(apr_1)
    Year1_monthly.append(may_1)
    Year1_monthly.append(jun_1)
    Year1_monthly.append(jul_1)
    Year1_monthly.append(aug_1)
    Year1_monthly.append(sep_1)
    Year1_monthly.append(oct_1)
    Year1_monthly.append(nov_1)
    Year1_monthly.append(dec_1)

    print(Year1_monthly)
    print(len(Year1_monthly))
    print('Load :', user_load)
    print(len(user_load))

# Apply load escalation annually and save every year

def annual_load_escalation(load_esc):
    annual_load = []
    #Applying escalation for 25 years
    for n in range(0,26):
        esc_load_n = user_load
        esc_load_n = esc_load_n * (1 + (n * load_esc))
        annual_load.append(esc_load_n)

    return annual_load

print('The annual load for year:', len(annual_load_escalation(load_esc)[1]))


# apply escalation to load and save in monthly for 25 years
def monthly_cum_load(n):
    monthly_cum_load = []
    esc_load = annual_load_escalation(load_esc)[n]
    jan_1 = 0
    feb_1 = 0
    mar_1 = 0
    apr_1 = 0
    may_1 = 0
    jun_1 = 0
    jul_1 = 0
    aug_1 = 0
    sep_1 = 0
    oct_1 = 0
    nov_1 = 0
    dec_1 = 0

# load cummulative for monthly
    for i in range(0, 8760):
        if i in range(0, 744):
            jan_1 = jan_1 + esc_load[i]
            monthly_cum_load.append(jan_1)
        elif i in range(745, 1416):
            feb_1 = feb_1 + esc_load[i]
            monthly_cum_load.append(feb_1)
        elif i in range(1417, 2160):
            mar_1 = mar_1 + esc_load[i]
            monthly_cum_load.append(mar_1)
        elif i in range(2161, 2880):
            apr_1 = apr_1 + esc_load[i]
            monthly_cum_load.append(apr_1)
        elif i in range(2881, 3624):
            may_1 = may_1 + esc_load[i]
            monthly_cum_load.append(may_1)
        elif i in range(3625, 4344):
            jun_1 = jun_1 + esc_load[i]
            monthly_cum_load.append(jun_1)
        elif i in range(4345, 5088):
            jul_1 = jul_1 + esc_load[i]
            monthly_cum_load.append(jul_1)
        elif i in range(5089, 5832):
            aug_1 = aug_1 + esc_load[i]
            monthly_cum_load.append(aug_1)
        elif i in range(5832, 6552):
            sep_1 = sep_1 + esc_load[i]
            monthly_cum_load.append(sep_1)
        elif i in range(6553, 7296):
            oct_1 = oct_1 + esc_load[i]
            monthly_cum_load.append(oct_1)
        elif i in range(7297, 8016):
            nov_1 = nov_1 + esc_load[i]
            monthly_cum_load.append(nov_1)
        else:
            dec_1 = dec_1 + esc_load[i]
            monthly_cum_load.append(dec_1)

    # print('The cumulative monthly load is:', len(monthly_cum_load))

    return monthly_cum_load

# print('The monthly escalated load:',monthly_cum_load(1))
# print('The total load for Year 2 :', sum(load_escalation(load_esc)[2]))

# selection of appropriate slab with respect to the monthly avg consumption
def slab_selection(n):
    avg_monthly = sum(annual_load_escalation(load_esc)[n]) / 12

    for s in range(1,slab_id_m+1 ):
        # print('Test')
        EC_t_q = pd.read_sql_query("select slab_id, period, tier, min, maximum, energy_charge from slabs_mapping where state_id =" + str(
                            state_id) + " and tarriff_type_id = " + str(
                            tariff_id) + " and voltage_type_id = " + str(
                            voltage_id) + " and metering_type_id = " + str(metering_id) + " and slab_id = " + str(
                            s), conn)
        EC_t = EC_t_q
        max_slab: float = EC_t['maximum'].max(0)
        # print('THe largest max is:', max_slab)
        min_slab: float = EC_t['min'].min(0)
        # print('The smallest min is:', min_slab)

        if float(max_slab) > float(avg_monthly) and float(avg_monthly) > float(min_slab):
            EC_matrix = EC_t
            break
           # print('The EC slab is:', EC_matrix)

    return EC_matrix

# print('The selected EC table is:', slab_selection(0))

# apply escalation to the fixed charge calculation
def fixed_charge_esc(n):
    FC = fixed_charge_h

    for j in range(n, n+1):
        FC_esc = FC * (1 + (j * cost_esc))

    # print('The escalated fixed charge is:', FC_esc)

    return FC_esc
# print('The fixed charge is:', fixed_charge_esc(0))

# apply escalation to cost and save in monthly for 25 years
def cost_escalation(n):
    EC_cost_esc = pd.DataFrame(slab_selection(n))
    temp_EC = EC_cost_esc['energy_charge']
    # load_esc_years = [5, 10, 15, 20, 25]

    for j in range(n, n+1):
        # if j in load_esc_years:
        # print('The table is:', EC_cost_esc)
        temp_EC = EC_cost_esc['energy_charge']
        EC_cost_esc['energy_charge'] = temp_EC * (1 + (j * cost_esc))
    slab_id_a = EC_cost_esc['slab_id'].max()


    return EC_cost_esc, slab_id_a

# print('The applicable energy charge is:',(cost_escalation(1)[0]))
# print('The applicable slab is:',(cost_escalation(1)[1]))
# print('The cost escalation for Year 2 is :' + str(cost_escalation(cost_energy, cost_esc)[1]))

# sorting 25 year monthly bill from Elec_bill_wo_system
# Build a TOD matrix of each month based on applicability and periods
# Check if TOD is applicable or not

# tou_select_q = pd.read_sql_query("select applicability_periods from tou_applicability_periods where tariff_type_id=" + str(tariff_id) + " and voltage_type_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
# tou_select = tou_select_q
# print(tou_select)
# Matrix for applying TOU
if tou_select == 0:
# Build matrix with 8760 points for normal period (i.e. = 1)
    tou_matrix = [1] * 8760
    # print(len(tou_matrix))
elif tou_select == 1 or tou_select == 2:
    tou_matrix = [1] * 24
    # print(tou_matrix)
#identify max rows applicable for the inputs
    tou_period_q = pd.read_sql_query("select tou_period from tou_period_table where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
    tou_period = tou_period_q
    # print('check:', tou_period)
    tou_from_q = pd.read_sql_query("select tou_from_hr from tou_period_table where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
    tou_from:float = tou_from_q
    # print('check from:', tou_from)
    tou_to_q = pd.read_sql_query("select tou_to_hr from tou_period_table where tariff_id=" + str(tariff_id) + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
    tou_to:float = tou_to_q
    # print('check from:', tou_to)
    tou_from_to = []
    tou_from_to = tou_period
    df = pd.DataFrame(tou_from_to)
    df['tou_from_hr'] = tou_from
    df['tou_to_hr'] = tou_to
    # print('Period:', df)
    # print(tou_matrix[0:10])
    # print(df['tou_from_hr'][0])
    # print(len(tou_period))
# Replacing the tou matrix with applicable peak (2) and off peak (3) periods
    # initialising range in tou matrix
    for u in range(len(tou_period)):
        if df['tou_period'][u] == 2:
            r1, r2 = df['tou_from_hr'][u], df['tou_to_hr'][u]
            r = 2
        elif df['tou_period'][u] == 3:
            r1, r2 = df['tou_from_hr'][u], df['tou_to_hr'][u]
            r = 3
        tou_matrix[r1:r2] = [r] * (r2 - r1)
    # print('TOU Matrix:',tou_matrix)

    tou_matrix = tou_matrix * 365
    # print(len(tou_matrix))

# Identify the TOU table
# Peak TOU table
TOU_p_q = pd.read_sql_query("select period, tier, min, maximum, energy_charge from tou_period_energy_charge where tariff_type_id=" + str(tariff_id) + " and voltage_type_id=" + str(voltage_id) + " and state_id=" + str(state_id) + " and period=" + str(2), conn)
TOU_p = TOU_p_q

# Off peak TOU table
TOU_op_q = pd.read_sql_query("select period, tier, min, maximum, energy_charge from tou_period_energy_charge where tariff_type_id=" + str(tariff_id) + " and voltage_type_id=" + str(voltage_id) + " and state_id=" + str(state_id) + " and period=" + str(3), conn)
TOU_op = TOU_op_q

# print('The peak TOU table is:', TOU_p)
# print('The off peak TOU table is:', TOU_op)

# Escalate the TOU peak charges
def cost_escalation_p(n):
    EC_p_cost_esc = pd.DataFrame(TOU_p)
    temp_EC_p = EC_p_cost_esc['energy_charge']
    # load_esc_years = [5, 10, 15, 20, 25]

    for j in range(n, n+1):
        # if j in load_esc_years:
            # print('The table is:', EC_cost_esc)
        temp_EC_p = EC_p_cost_esc['energy_charge']
        EC_p_cost_esc['energy_charge'] = temp_EC_p * (1 + (j * cost_esc))
    # slab_id_a = EC_p_cost_esc['slab_id'].max()
    return EC_p_cost_esc
# print('The applicable peak energy charge is:',(cost_escalation_p(1)))

# Escalate the TOU off peak charges
def cost_escalation_op(n):
    EC_op_cost_esc = pd.DataFrame(TOU_op)
    temp_EC_op = EC_op_cost_esc['energy_charge']
    # load_esc_years = [5, 10, 15, 20, 25]

    for j in range(n, n+1):
        # if j in load_esc_years:
            # print('The table is:', EC_cost_esc)
        temp_EC_op = EC_op_cost_esc['energy_charge']
        EC_op_cost_esc['energy_charge'] = temp_EC_op * (1 + (j * cost_esc))
    # slab_id_a = EC_p_cost_esc['slab_id'].max()
    return EC_op_cost_esc
# print('The applicable off peak energy charge is:',(cost_escalation_op(1)))


# Selection of applicable charges
# Build a cost matrix with 8760 points
def charge_selection_n(n):
    TOU_matrix = tou_matrix
    cum_monthly_load = monthly_cum_load(n)
    # print('THe cummulative load is:', cum_monthly_load)
    app_EC = cost_escalation(n)[0]
    app_EC_p = cost_escalation_p(n)
    app_EC_op = cost_escalation_op(n)
    # print('The applied EC is:',app_EC)
    # slab_id_app = cost_escalation(n)[1]
    tier_app = app_EC['tier'].max()
    cost_matrix = []

    for i in range(0, 8760):
        if TOU_matrix[i] == 1:
            for t in range(1, tier_app + 1):
                if (app_EC['maximum'][t - 1]) >= (cum_monthly_load[i]) and (cum_monthly_load[i]) >= (app_EC['min'][t - 1]):
                    break
            cost_energy = app_EC['energy_charge'][t - 1]

        elif TOU_matrix[i] == 2:
            for t in range(1, tier_app + 1):
                    if (app_EC_p['maximum'][t - 1]) >= (cum_monthly_load[i]) and (cum_monthly_load[i]) >= (app_EC_p['min'][t - 1]):
                        break
            cost_energy = app_EC_p['energy_charge'][t - 1]

        elif TOU_matrix[i] == 3:
            for t in range(1, tier_app + 1):
                    if (app_EC_op['maximum'][t - 1]) >= (cum_monthly_load[i]) and (cum_monthly_load[i]) >= (app_EC_op['min'][t - 1]):
                        break
            cost_energy = app_EC_op['energy_charge'][t - 1]

        cost_matrix.append(cost_energy)


    # print('The cost matrix is:', (cost_matrix))

    return cost_matrix

# print('The energy cost matrix is:',charge_selection_n(1))


#Calculate the bill without system

def Elec_bill_wo_system(n):
    Elec_bill_wo_system = 0
    bill_load = annual_load_escalation(load_esc)[n]
    bill_cost = charge_selection_n(n)
    fixed_charge_h_n = fixed_charge_esc(n)
    bill_annual = 0

    for i in range(0, 8760):
        bill_annual = bill_annual + (fixed_charge_h_n + (bill_load[i] * bill_cost[i]))

    Elec_bill_wo_system = bill_annual

    # print('The annual bill w/o system is:', Elec_bill_wo_system)

    return Elec_bill_wo_system

print('The annual bill is:', Elec_bill_wo_system(0))


# power balance for the selected size
def power_balance(x1, g):
    # print(x1[0])
    # print(type(x1))
    # print(n)
    rem_load = np.zeros(n)
    exsolar = np.zeros(n)
    sol_cap = x1[0]
    # print(type(sol_cap))
    # print(sol_cap)
    # print(type(solarp))
    # print((solarp['load'][10]))


    # print('The solar capacity is:',sol_cap)
    solar_power = (solarp['load'] * sol_cap) * (1 - (g * der_deg))
    # print((solar_power[10]))

    # print('The solar power generated:', sum(solar_power))
    # print(len(solar_power))
    bat_cap = x1[1]
    # print('The battery capacity is:',bat_cap)
    bat_inv = 0.835 * sol_cap
    ch_dis_available = np.zeros(n)
    soc = np.zeros(n)

    battery_power = np.zeros(n)
    gridp = np.zeros(n)
    gridpower = np.zeros(n)
    excessder = np.zeros(n)
    load = annual_load_escalation(load_esc)[g]
    # print(len(user_load))
    # print('The load is :', load)

    #Power balance for Gross metering
    for i in range(0, 8760):
        # print('Load', user_load[i])
        # print('Solar', solar_power[i][0])
        # print(type(solar_power[i][0]))

    #Identifying the available remaining load & excess solar
        # if metering_type == "Gross Metering":
        #     rem_load[i] = load[i]
        #
        #     exsolar[i] = (solar_power[i])
        #     # print('Remaining load:', sum(rem_load))
        #     # print('Excess solar:', sum(exsolar))

        if (load[i] - (solar_power[i])) > 0:
            rem_load[i] = load[i] - solar_power[i]
        else:
            exsolar[i] = - load[i] + solar_power[i]

    # calculating charging/discharging power available
    if solar & battery:
        for i in range(n):
            if batstatus[0][i] == -1:
                if exsolar[i] < bat_inv:
                    ch_dis_available[i] = exsolar[i]
                else:
                    ch_dis_available[i] = bat_inv

            elif batstatus[0][i] == 1:
                if rem_load[i] < bat_inv:
                    ch_dis_available[i] = rem_load[i]
                else:
                    ch_dis_available[i] = bat_inv

        # calculating the state of charge based on min and max limit of soc
        socbatmax = socmax * bat_cap
        # print('SOC max', socbatmax)
        socbatmin = socmin * bat_cap
        # print('SOC min', socbatmin)
        soc[0] = socin * bat_cap
        # print('SOC initial', soc[0])

        for i in range(n - 1):
            k = i + 1
            if batstatus[0][k] == -1:
                if (ch_dis_available[k] + soc[k - 1]) < socbatmax:
                    soc[k] = ch_dis_available[k] + soc[k - 1]
                else:
                    soc[k] = socbatmax
            if batstatus[0][k] == 1:
                if (soc[k - 1] - ch_dis_available[k]) > socbatmin:
                    soc[k] = soc[k - 1] - ch_dis_available[k]
                else:
                    soc[k] = socbatmin
            elif batstatus[0][k] == 0:
                soc[k] = soc[k - 1]

            battery_power[k] = soc[k - 1] - soc[k]  # calculating actual battery discharging and charging power
            # print('Battery power:', (battery_power))
        for i in range(0, 8760):
            if metering_type == "Gross Metering":
                gridp[i] = load[i]
                gridpower[i] = gridp[i]
                excessder[i] = (solar_power[i]) - battery_power[i]

            else :
                gridp[i] = load[i] - (solar_power[i]) - battery_power[i]  # allocating export and grid after battery
                if gridp[i] > 0:
                    gridpower[i] = gridp[i]
                    excessder[i] = 0
                elif gridp[i] < 0:
                    gridpower[i] = 0
                    excessder[i] = -gridp[i] # what happens if battery has remaining energy? Does it go to excess DER?
            # print(soc)
            # print(battery_power)
    else:
        if metering_type == "Gross Metering":
            gridpower = load
            excessder = solar_power
        else:
            gridpower = rem_load
            excessder = exsolar


    df1 = pd.DataFrame()
    df1['load'] = pd.DataFrame(load)
    df1['solar'] = solar_power
    df1['grid'] = gridpower
    df1['battery'] = battery_power
    df1['excess'] = excessder
    # df1['SOC'] = soc[k]

    sum_solar = sum(solar_power)
    sum_battery = sum(battery_power)
    sum_grid = sum(gridpower)
    sum_export = sum(excessder)
    if metering_type == "Gross Metering":
        sum_solar_load = 0
    else:
        sum_solar_load = sum_solar - sum_battery - sum_export


    # print('daily solar :', sum_solar/365)
    # print('Export:', sum(excessder))
    return battery_power, gridpower, excessder, solar_power, load, df1, sum_battery, sum_grid, sum_export, sum_solar, sum_solar_load


# print('The power balance details are:', (power_balance(x1, 0)[5]))

# print('The total solar generation is:', power_balance(x1, 0)[9])
# print('The total solar contribution is:', power_balance(x1, 0)[10])
# print('The total battery contribution is:', power_balance(x1, 0)[6])
# print('The total grid contribution is:', power_balance(x1, 0)[7])
# print('The total export contribution is:', power_balance(x1, 0)[8])


# apply escalation to grid and export also save in monthly for 25 years
def grid_monthly_cum_load(x1, n):
    g_monthly_cum_load = []
    g_esc_load = power_balance(x1, n)[1]
    e_monthly_cum_load = []
    e_esc_load = power_balance(x1, n)[2]
# variables to store grid values
    g_jan_1 = 0
    e_jan_1 = 0
    g_feb_1 = 0
    e_feb_1 = 0
    g_mar_1 = 0
    e_mar_1 = 0
    g_apr_1 = 0
    e_apr_1 = 0
    g_may_1 = 0
    e_may_1 = 0
    g_jun_1 = 0
    e_jun_1 = 0
    g_jul_1 = 0
    e_jul_1 = 0
    g_aug_1 = 0
    e_aug_1 = 0
    g_sep_1 = 0
    e_sep_1 = 0
    g_oct_1 = 0
    e_oct_1 = 0
    g_nov_1 = 0
    e_nov_1 = 0
    g_dec_1 = 0
    e_dec_1 = 0


# load cummulative for monthly
    for i in range(0, 8760):
        if i in range(0, 744):
            g_jan_1 = g_jan_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_jan_1)
            e_jan_1 = e_jan_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_jan_1)
        elif i in range(745, 1416):
            g_feb_1 = g_feb_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_feb_1)
            e_feb_1 = e_feb_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_feb_1)
        elif i in range(1417, 2160):
            g_mar_1 = g_mar_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_mar_1)
            e_mar_1 = e_mar_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_mar_1)
        elif i in range(2161, 2880):
            g_apr_1 = g_apr_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_apr_1)
            e_apr_1 = e_apr_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_apr_1)
        elif i in range(2881, 3624):
            g_may_1 = g_may_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_may_1)
            e_may_1 = e_may_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_may_1)
        elif i in range(3625, 4344):
            g_jun_1 = g_jun_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_jun_1)
            e_jun_1 = e_jun_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_jun_1)
        elif i in range(4345, 5088):
            g_jul_1 = g_jul_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_jul_1)
            e_jul_1 = e_jul_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_jul_1)
        elif i in range(5089, 5832):
            g_aug_1 = g_aug_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_aug_1)
            e_aug_1 = e_aug_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_aug_1)
        elif i in range(5832, 6552):
            g_sep_1 = g_sep_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_sep_1)
            e_sep_1 = e_sep_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_sep_1)
        elif i in range(6553, 7296):
            g_oct_1 = g_oct_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_oct_1)
            e_oct_1 = e_oct_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_oct_1)
        elif i in range(7297, 8016):
            g_nov_1 = g_nov_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_nov_1)
            e_nov_1 = e_nov_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_nov_1)
        else:
            g_dec_1 = g_dec_1 + g_esc_load[i]
            g_monthly_cum_load.append(g_dec_1)
            e_dec_1 = e_dec_1 + e_esc_load[i]
            e_monthly_cum_load.append(e_dec_1)

    # print('The cumulative monthly load is:', len(g_monthly_cum_load))

    return g_monthly_cum_load, e_monthly_cum_load
# print('The cumulative export is:', (grid_monthly_cum_load(0)[1]))

# slab selection after system based on avg monthly consumption
def g_slab_selection(n):
    g_avg_monthly = sum(annual_load_escalation(load_esc)[n]) / 12

    for s in range(1,slab_id_m+1 ):
        # print('Test')
        g_EC_t_q = pd.read_sql_query("select slab_id, period, tier, min, maximum, energy_charge from slabs_mapping where state_id =" + str(
                            state_id) + " and tarriff_type_id = " + str(
                            tariff_id) + " and voltage_type_id = " + str(
                            voltage_id) + " and metering_type_id = " + str(metering_id) + " and slab_id = " + str(
                            s), conn)
        g_EC_t = g_EC_t_q
        g_max_slab: float = g_EC_t['maximum'].max(0)
        # print('THe largest max is:', max_slab)
        g_min_slab: float = g_EC_t['min'].min(0)
        # print('The smallest min is:', min_slab)

        if float(g_max_slab) > float(g_avg_monthly) and float(g_avg_monthly) > float(g_min_slab):
            g_EC_matrix = g_EC_t
            break
           # print('The EC slab is:', EC_matrix)

    return g_EC_matrix

# print('The selected EC table is:', slab_selection(0))

# apply escalation to the energy charge for 25 years
def g_cost_escalation(n):
    g_EC_cost_esc = pd.DataFrame(g_slab_selection(n))
    g_temp_EC = g_EC_cost_esc['energy_charge']
    # load_esc_years = [5, 10, 15, 20, 25]

    for j in range(n, n+1):
        # if j in load_esc_years:
        # print('The table is:', EC_cost_esc)
        g_temp_EC = g_EC_cost_esc['energy_charge']
        g_EC_cost_esc['energy_charge'] = g_temp_EC * (1 + (j * cost_esc))
    g_slab_id_a = g_EC_cost_esc['slab_id'].max()


    return g_EC_cost_esc, g_slab_id_a

# print('The applicable energy charge is:',(g_cost_escalation(5)[0]))

# Selection of applicable charges
# Build a cost matrix with 8760 points
def g_charge_selection_n(x1, n):
    g_TOU_matrix = tou_matrix
    g_cum_monthly_load = grid_monthly_cum_load(x1, n)[0]
    # print('THe cummulative load is:', cum_monthly_load)
    g_app_EC = g_cost_escalation(n)[0]
    g_app_EC_p = cost_escalation_p(n)
    g_app_EC_op = cost_escalation_op(n)
    # print('The applied EC is:',app_EC)
    # slab_id_app = cost_escalation(n)[1]
    g_tier_app = g_app_EC['tier'].max()
    g_cost_matrix = []

    for i in range(0, 8760):
        if g_TOU_matrix[i] == 1:
            for t in range(1, g_tier_app + 1):
                if (g_app_EC['maximum'][t - 1]) >= (g_cum_monthly_load[i]) and (g_cum_monthly_load[i]) >= (g_app_EC['min'][t - 1]):
                    break
            g_cost_energy = g_app_EC['energy_charge'][t - 1]

        elif g_TOU_matrix[i] == 2:
            for t in range(1, g_tier_app + 1):
                    if (g_app_EC_p['maximum'][t - 1]) >= (g_cum_monthly_load[i]) and (g_cum_monthly_load[i]) >= (g_app_EC_p['min'][t - 1]):
                        break
            g_cost_energy = g_app_EC_p['energy_charge'][t - 1]

        elif g_TOU_matrix[i] == 3:
            for t in range(1, g_tier_app + 1):
                    if (g_app_EC_op['maximum'][t - 1]) >= (g_cum_monthly_load[i]) and (g_cum_monthly_load[i]) >= (g_app_EC_op['min'][t - 1]):
                        break
            g_cost_energy = g_app_EC_op['energy_charge'][t - 1]

        g_cost_matrix.append(g_cost_energy)


    # print('The cost matrix is:', (g_cost_matrix))

    return g_cost_matrix

# print('The energy cost matrix is:',charge_selection_n(0))

# Network charges and escalation
def network_charge_esc(x1, n):
    NC_temp = float(network_charge_fetch(x1[0])[0])
    # load_esc_years = [5, 10, 15, 20, 25]

    for n in range(n, n+1):
        # if n in load_esc_years:
        NC_temp = NC_temp * (1 + (n * cost_esc))

    NC_cost = NC_temp
    # print('The network charge is:', str(NC_cost))

    return NC_cost

# print('The network charge is:', network_charge_esc(0))

# Sorting normal, peak and off peak for each month

def TOU_sort(x1, n):
    grid_t = power_balance(x1, n)[1]
    g_units_n = [0]*8760
    g_units_p = [0]*8760
    g_units_op = [0]*8760
    g_unit_n_s = 0
    g_unit_p_s = 0
    g_unit_op_s = 0
    g_units_n_s = []
    g_units_p_s = []
    g_units_op_s = []
    end_of_month = [0, 744, 1416, 2160, 2880, 3624, 4344, 5088, 5832, 6552, 7296, 8016, 8760]
    # print('The month index is:',end_of_month[2])
    # print('The length are:', len(end_of_month))
    for i in range(0,8760):
        if tou_matrix[i] == 1:
            g_units_n[i] = grid_t[i]
        elif tou_matrix[i] == 2:
            g_units_p[i] = grid_t[i]
        elif tou_matrix[i] == 3:
            g_units_op[i] = grid_t[i]
    #
    # g_unit_n_s = sum(g_units_n[end_of_month[0]:end_of_month[1]])
    # print('The sum is', g_unit_n_s)

    #sum of normal peak and off peak for each month
    for i in range(0,12):
        g_unit_n_s = sum(g_units_n[end_of_month[i]:end_of_month[i+1]])
        g_units_n_s.append(g_unit_n_s)
        g_unit_p_s = sum(g_units_p[end_of_month[i]:end_of_month[i + 1]])
        g_units_p_s.append(g_unit_p_s)
        g_unit_op_s = sum(g_units_op[end_of_month[i]:end_of_month[i + 1]])
        g_units_op_s.append(g_unit_op_s)

    # print('The TOU normal is:', (g_units_n_s))
    # print('The TOU peak is:', (g_units_p_s))
    # print('The TOU off peak is:', (g_units_op_s))

    return g_units_n_s, g_units_p_s, g_units_op_s
# print('The TOU split is:', TOU_sort(0))

# Sorting grid and export units for Net Metering billing
def cum_export_unit(x1, n):
    g_units = []
    e_units = []
    d_units = []
    d_units_n = []
    e_units_t = 0
    E_units = []
    grid_u = grid_monthly_cum_load(x1, n)[0]
    # print(len(grid_u))
    export_u = grid_monthly_cum_load(x1, n)[1]
    end_of_month = [743, 1415, 2159, 2879, 3623, 4343, 5087, 5831, 6551, 7295, 8015, 8759]
    df_net = pd.DataFrame()

    for i in range(0,8760):
        d_unit = grid_u[i] - export_u[i]
        if i in end_of_month:
            g_units.append(grid_u[i])
            e_units.append(export_u[i])
            d_units.append(d_unit)
    # print('The length of matrix is:',len(d_units))
    # for i in range(0,12):
    #     d_unit = g_units[i] - e_units[i]
    #     d_units.append(d_unit)

    d_units_n.append(d_units[0])
    E_units.append(e_units[0])
    for i in range(1,12):
        if d_units_n[i-1] >= 0:
            E_units.append(e_units[i])
            d_unit_n = g_units[i] - E_units[i]
            d_units_n.append(d_unit_n)
        else:
            e_units_t = e_units[i] - d_units_n[i-1]
            E_units.append(e_units_t)
            d_unit_n = g_units[i] - E_units[i]
            d_units_n.append(d_unit_n)
    e_units_n = [ i * -1 for i in e_units] # applying a negative to the export

    # df_net['Grid total'] = pd.DataFrame(g_units)
    df_net['Grid normal'] = TOU_sort(x1, n)[0]
    df_net['Grid peak'] = TOU_sort(x1, n)[1]
    df_net['Grid off peak'] = TOU_sort(x1, n)[2]
    # df_net['Excess units'] = e_units_n
    # df_net['Difference of units'] = d_units
    # df_net['Updated Excess'] = E_units
    df_net['Updated Difference'] = d_units_n

    # print('The monthly cumulative units:', (df_net))
    # print('The excess units', (e_units))
    # print('The diff units is:', (d_units))
    # print('The updated excess units is:', (E_units))
    # print('The updted diff units:', (d_units_n))

    return df_net
# print('The grid unit are:',cum_export_unit(0))

# Sorting billable units for Net metering
def NM_bill(x1, n):
    diff_up = cum_export_unit(x1, n)['Updated Difference']
    peak_u = cum_export_unit(x1, n)['Grid peak']
    offpeak_u = cum_export_unit(x1, n)['Grid off peak']
    normal_u = cum_export_unit(x1, n)['Grid normal']
    bill_units_p = [0]*12
    # print(len(bill_units_p))
    bill_units_op = [0]*12
    # print(len(bill_units_op))
    bill_units_n = [0]*12
    # print(len(diff_up))
    bill_unit_p = 0
    bill_unit_n = 0
    bill_unit_op = 0
    NM_df = pd.DataFrame()

    for i in range(0,12):
        if diff_up[i] > 0:
            bill_unit_p = diff_up[i] - peak_u[i]
            if bill_unit_p > 0:
                bill_units_p[i] = peak_u[i]
                bill_unit_n = bill_unit_p - normal_u[i]
            else:
                bill_units_p[i] = diff_up[i]


            if bill_unit_n > 0:
                bill_units_n[i] = normal_u[i]
                bill_unit_op = bill_unit_n - offpeak_u[i]
            else:
                bill_units_n[i] = bill_unit_p


            if bill_unit_op > 0:
                bill_units_op[i] = offpeak_u[i]
            else:
                bill_units_op[i] = bill_unit_n

        else:
            bill_units_p[i] = 0
            bill_units_n[i] = 0
            bill_units_op[i] = 0



    NM_df['Peak'] = bill_units_p
    NM_df['Normal'] = bill_units_n
    NM_df['Offpeak'] = bill_units_op
    NM_df.insert(0, 'Diff', diff_up)# Indexing it to first column

    NM_df = NM_df.where(NM_df > 0, 0)# replacing negative values to zero

    # NM_df.insert(0, 'Diff', diff_up)
    # NM_df['Diff'] = diff_up

    # print('The updated diff is:', NM_df)

    return NM_df
# print('The NM bill is:', NM_bill(0))
# print(type(NM_bill(0)))

# billable units are stored here
def nm_monthly_cal(x1, n):

    total_bill: float = NM_bill(x1, n)['Diff']
    # print('The total bill is :', total_bill)
    return total_bill

# print('The monthly bill is:', nm_monthly_cal(0))

# slab selection based on billable units
def nm_slab_selection(x1, n):
    nm_avg_monthly: float = nm_monthly_cal(x1, n)
    # print(type(nm_avg_monthly))
    nm_slab = [0]*12

    for i in range(len(nm_avg_monthly)):

        for s in range(1,slab_id_m+1 ):
            # print('Test')
            nm_EC_t_q = pd.read_sql_query("select slab_id, period, tier, min, maximum, energy_charge from slabs_mapping where state_id =" + str(
                                state_id) + " and tarriff_type_id = " + str(
                                tariff_id) + " and voltage_type_id = " + str(
                                voltage_id) + " and metering_type_id = " + str(metering_id) + " and slab_id = " + str(
                                s), conn)
            nm_EC_t = nm_EC_t_q
            nm_max_slab: float = nm_EC_t['maximum'].max(0)
            # print('THe largest max is:', type(nm_max_slab))
            nm_min_slab: float = nm_EC_t['min'].min(0)
            # print('The smallest min is:', type(nm_min_slab))


            if (nm_max_slab) > float (nm_avg_monthly[i]) and float (nm_avg_monthly[i]) > (nm_min_slab):
                nm_EC_matrix = nm_EC_t
                nm_slab[i] = s
                break
        # print('The EC slab is:', nm_slab)

    return nm_EC_matrix, nm_slab

# print('The selected EC table is:', nm_slab_selection(0)[1])

# NM bill calculation
def NM_bill_cal(x1, n):
    unit_matrix = NM_bill(x1, n)
    EC_n = slab_id_t
    EC_p = TOU_p
    EC_op = TOU_op
    slab_matrix = nm_slab_selection(x1, n)[1]
    unit_matrix = unit_matrix.drop('Diff', axis=1)
    Total_bill = pd.DataFrame()
    Total_bill['Bill'] = [0]*12
# Check each column of row to calculate peak, normal and off peak charges
    for i in range (0,12):
        # Total_bill['Bill'][i] = 5
        bill_t_n: float = 0
        bill_t_p: float = 0
        bill_t_op: float = 0
        n_total = 0
        p_total = 0
        op_total = 0
        # print('The iteration is:', i)
        unit_matrix_t = unit_matrix.iloc[i:i+1]
# peak units calculation
        if unit_matrix_t['Peak'].all() > 0:
            EC_p_t = EC_p
            m_tier = EC_p_t['tier'].max()
            if m_tier > 1:
                EC_p_t = EC_p.loc[EC_p['period'] == 2]
                EC_p_t = EC_p_t.reset_index()

            m_tier_p = EC_p_t['tier'].max()
            # print('The max tier is:', m_tier)
            # print('The first value of max is',EC_n_t['maximum'])
            units_p: float = unit_matrix_t['Peak']

            for j in range(1,m_tier_p+1):
                # print('The units are:', bill_t_n)
                # print('The first element is:',float(units_n))
                if float(units_p) >= float(EC_p_t['maximum'][j-1]):
                    bill_t_p = bill_t_p + ((EC_p_t['maximum'][j - 1] - EC_p_t['min'][j - 1]) * EC_p_t['energy_charge'][j - 1])
                    units_p = units_p - (EC_p_t['maximum'][j - 1] - EC_p_t['min'][j - 1])
                    units_p = units_p
                    # print('the unit is:', units_n)
                else:
                    bill_t_p = bill_t_p + (units_p * EC_p_t['energy_charge'][j-1])
                    units_p = units_p - units_p

            # print('The peak bill amount is:', bill_t_p)
                # print('The temp units are:', units_n)
            p_total = bill_t_p

# normal units calculation
        if unit_matrix_t['Normal'].all() > 0:
            EC_n_t = EC_n
            m_tier = EC_n_t['tier'].max()
            if m_tier > 1:
                EC_n_t = EC_n.loc[EC_n['slab_id'] == slab_matrix[i]]
                EC_n_t = EC_n_t.reset_index()

            m_tier_n = EC_n_t['tier'].max()
            # print('The max tier is:', m_tier)
            # print('The first value of max is',EC_n_t['maximum'])
            units_n: float = unit_matrix_t['Normal']

            for j in range(1,m_tier_n+1):
                # print('The units are:', bill_t_n)
                # print('The first element is:',float(units_n))
                if float(units_n) >= float(EC_n_t['maximum'][j-1]):
                    bill_t_n = bill_t_n + ((EC_n_t['maximum'][j - 1] - EC_n_t['min'][j - 1]) * EC_n_t['energy_charge'][j - 1])
                    units_n = units_n - (EC_n_t['maximum'][j - 1] - EC_n_t['min'][j - 1])
                    units_n = units_n
                    # print('the unit is:', units_n)
                else:
                    bill_t_n = bill_t_n + (units_n * EC_n_t['energy_charge'][j-1])
                    units_n = units_n - units_n

            # print('The bill amount is:', bill_t_n)
                # print('The temp units are:', units_n)
            n_total = bill_t_n
# offpeak units calculation
        if unit_matrix_t['Offpeak'].all() > 0:
            EC_op_t = EC_op
            m_tier = EC_op_t['tier'].max()
            if m_tier > 1:
                EC_op_t = EC_op.loc[EC_op['period'] == 3]
                EC_op_t = EC_op_t.reset_index()

            m_tier_op = EC_op_t['tier'].max()
            # print('The max tier is:', m_tier)
            # print('The first value of max is',EC_n_t['maximum'])
            units_op: float = unit_matrix_t['Offpeak']

            for j in range(1, m_tier_op + 1):
                # print('The units are:', bill_t_n)
                # print('The first element is:',float(units_n))
                if float(units_op) >= float(EC_op_t['maximum'][j - 1]):
                    bill_t_op = bill_t_op + (
                                (EC_op_t['maximum'][j - 1] - EC_op_t['min'][j - 1]) * EC_op_t['energy_charge'][j - 1])
                    units_op = units_op - (EC_op_t['maximum'][j - 1] - EC_op_t['min'][j - 1])
                    units_op = units_op
                    # print('the unit is:', units_n)
                else:
                    bill_t_op = bill_t_op + (units_op * EC_op_t['energy_charge'][j - 1])
                    units_op = units_op - units_op

                # print('The offpeak bill amount is:', bill_t_op)
                # print('The temp units are:', units_n)
                op_total = bill_t_op

        Total_bill['Bill'][i] = (n_total + p_total + op_total)
    # print('The total bill is :', sum(Total_bill['Bill']))

    NM_bill_t = sum(Total_bill['Bill'])

        # print((unit_matrix_t['Normal']))

    # print('The unit matrix :', unit_matrix)
    # print('The EC normal:', EC_n)
    # print('The EC peak:', EC_p)
    # print('The EC off peak:', EC_op)
    # print('The slab matrix :', slab_matrix)

    return NM_bill_t
# print('The NM bill is:', NM_bill_cal(0))

# bill with system
def Elec_bill_w_sys(x1, n):
# bill with system for Gross metering
    if metering_type == "Gross Metering":
        cons_grid = power_balance(x1, n)[1]
        # print('The length of array is', len(cons_grid))
        cons_solar = power_balance(x1, n)[3]
        cost_grid = charge_selection_n(n)
        cost_solar = network_charge_fetch(x1[0])[1]
        # print('The feed in tariff is:', cost_solar)
        fixed_charge_h_n = fixed_charge_esc(n)
        network_charge = network_charge_esc(x1, n)[0]
        # print('The NC is:', network_charge)
        bill_w_sys = 0
        g_bill_annual = 0

        for g in range(0, 8760):
            g_bill_annual = g_bill_annual + (fixed_charge_h_n + (((cons_grid[g] * cost_grid[g]) - (cons_solar[g] * cost_solar) + (cons_solar[g] * network_charge))))

        bill_w_sys = g_bill_annual
        # print('The bill without system (Gross metering) is:', sum(bill_w_sys))

    elif metering_type == "Net Feed In":
        cons_grid = power_balance(x1, n)[1]
        cons_solar = power_balance(x1, n)[3]
        cons_export = power_balance(x1, n)[2]
        cost_grid = charge_selection_n(n)
        cost_solar = network_charge_fetch(x1[0])[1]
        fixed_charge_h_n = fixed_charge_esc(n)
        network_charge = network_charge_esc(x1, n)
        # print('The NC is:', network_charge)
        bill_w_sys = 0
        nf_bill_annual = 0

        for g in range(0, 8760):
            nf_bill_annual = nf_bill_annual + (fixed_charge_h_n + (((cons_grid[g] * cost_grid[g]) - (cons_export[g] * cost_solar) + (cons_solar[g] * network_charge))))

        bill_w_sys = nf_bill_annual

    elif metering_type == "Net Metering":
        cons_solar = power_balance(x1, n)[3]
        fixed_charge_h_n = fixed_charge_esc(n)
        network_charge = network_charge_esc(x1, n)
        # print('The NC is:', network_charge)
        # bill_w_sys = 0
        # nm_bill_annual = 0
        NC_annual = 0

        FC_annual = fixed_charge_h_n * 8760
        EC_annual = NM_bill_cal(x1, n)

        for g in range(0, 8760):
            NC_annual = NC_annual + ((cons_solar[g] * network_charge))

        nm_bill_annual = (FC_annual + EC_annual + NC_annual)

        bill_w_sys = nm_bill_annual



    return bill_w_sys

# print(' The bill with system is:', Elec_bill_w_sys(x1, 0))

# Financial calculation for the selected system
def financial_calc(x1):
    # print('x1', x1)
    # print(type(x1))
    load = sum(power_balance(x1, 0)[4])
    batpr = power_balance(x1, 0)[6]
    gridprin = power_balance(x1, 0)[7]
    exgrid = power_balance(x1, 0)[8]
    solpower_s = power_balance(x1, 0)[10]
    solpower_t = power_balance(x1, 0)[9]
    # Fraction of load
    solar_f = solpower_s/load
    # print('The solar contribution is:',solar_f)
    battery_f = batpr/load
    # print('The battery contribution is:', battery_f)
    grid_f = gridprin/load
    # print('The grid contribution is:', grid_f)
    export_f = exgrid/solpower_t
    # print('The export contribution is:', export_f)
    # # amount invested calculation
    sol_cap = x1[0]
    bat_cap = x1[1]
    amount_invested = investmentcost_calculate(sol_cap, bat_cap)
    print('The total installation cost is:', amount_invested)
    rep_batterycost, rep_invertercost = replacement_cost(sol_cap, bat_cap)
    loan_principal_amount = amount_invested * pysam_debt_fraction / 100
    eq_amount = amount_invested * (1 - (pysam_debt_fraction / 100))
    # average monthly/yearly cash flow
    # for base year

    cCF_t = np.zeros(nyr)
    CF = np.zeros(nyr)
    Elec_bill_withoutDER = np.zeros(nyr)
    fc_yearly_charge = np.zeros(nyr)
    Elec_bill_withDER = np.zeros(nyr)
    total_savings = np.zeros(nyr)
    total_op_cost = np.zeros(nyr)
    total_cost = np.zeros(nyr)
    total_debt_yearly = np.zeros(nyr)
    NPV_to_Savings: float = 0

    # fc_yearly_charge[0] = 12 * fixed_charge
    Elec_bill_withoutDER[0] = Elec_bill_wo_system(0)
    # # Elec_bill_withDER[k] = (sum(gridpr) - sum(exgrid)) * cost_energy * (1 + (k * cost_esc))
    Elec_bill_withDER[0] = 0
    total_op_cost[0] = 0# no escalation for O&M, must include
    total_cost[0] = eq_amount
    total_savings[0] = 0
    CF[0] = (total_savings[0]) - total_cost[0]
    cCF = CF[0]
    cCF_t[0] = CF[0]

    # print('Elec_bill_withoutDER ', Elec_bill_withoutDER[0])
    # print('Elec_bill_withDER ', Elec_bill_withDER[0])
    # print('Solar Power : ', sum(solpower))
    # print('Grid Power : ', sum(gridprin))
    # print('Export Power : ', sum(exgrid))
    # print('Battery Power : ', sum(batpr))
    # print('cash flow[0] :', CF[0])
    # print('total savings[0] :', total_savings[0])
    # print('total op cost :', total_op_cost[0])
    # print('total cost :', total_cost[0])
    # print('eq amt :', eq_amount)

    emi = npf.pmt(loan_rate / 12, 12 * loan_period, -loan_principal_amount, 0)
    cum_cashflow: float = 0
    # load_k = user_load
    # gridpr = gridprin
    # load_esc_years = [5, 10, 15, 20, 25]
    # load_esc_years = [4, 9, 14, 19, 24]
    for i in range(1,26):
        # k = i + 1
        # if k in load_esc_years:
        #     load_k = annual_load_escalation(load_esc)[i]
        #     # gridpr = gridprin * (1 + k * load_esc)
        #     # solpower_k = solpower * (1 - k * der_deg)
        #     batpr = power_balance(x1, i)
        #     gridprin = power_balance(x1, i)
        #     exgrid = power_balance(x1, i)
        #     solpower = power_balance(x1, i)
        # fc_yearly_charge[k] = 12 * fixed_charge * (1 + k * cost_esc)
        Elec_bill_withoutDER[i] = Elec_bill_wo_system(i-1)
        # Elec_bill_withDER[k] = (sum(gridpr) - sum(exgrid)) * cost_energy * (1 + (k * cost_esc))
        # Elec_bill_withDER[i] = Elec_bill_w_sys(x1, i-1)
        total_op_cost[i] = (500 * sol_cap + 500 * bat_cap) * (1 + ((i-1) * cost_esc))
        # print('Total op cost year 2:', total_op_cost[2])
        total_savings[i] = (Elec_bill_withoutDER[i] - Elec_bill_w_sys(x1, i-1))
        # print('Total saving:', total_savings[1])

        if i in range(2,12):
            total_debt_yearly[i] = 12 * emi
        total_cost[i] = total_op_cost[i] + total_debt_yearly[i] + rep_invertercost[i] + rep_batterycost[i]
        # print('Inv, replacement :', sum(rep_invertercost))
        # print('Bat, replacement :', sum(rep_batterycost))
        CF[i] = total_savings[i] - total_cost[i]
        # print('Total saving:', total_savings[1])
        # print('Total cost:', total_cost[1])
        # print(CF)
        cCF = cCF + CF[i]
        cCF_t[i] = cCF
        # print('total_savings[0]',  total_savings[0])

        cum_cashflow = cum_cashflow + CF[i]
    # print('The cumulative CF is:', cCF_t)

    bau_npv: float = npf.npv(dis_factor, Elec_bill_withoutDER[1:n])
    bau_npv = (bau_npv + Elec_bill_withoutDER[0])

    dis_saving: float = npf.npv(dis_factor, total_savings[1:n])
    dis_saving = (dis_saving + total_savings[0])

    npv = npf.npv(dis_factor, CF[1:n])
    # print('NPV[1,24]:', npv)
    npv = npv + CF[0]
    # print('NPV:', npv)

    NPV_to_Savings = ((npv / bau_npv) * 100)
    NPV_to_Savings = NPV_to_Savings

    # print(npv)
    cum_cashflow = cum_cashflow + CF[0]
    # print('cumcashflow:', cum_cashflow)

    average_annualcashflow = cum_cashflow / nyr
    # print('average_annualcashflow:', average_annualcashflow)
    # print(Elec_bill_withDER)
    # print(Elec_bill_withoutDER)
    # print(total_cost)
    # print(CF)
    # print('Elec. bill with sys Yr1', Elec_bill_withDER[0])
    # print('Elec. bill w/o sys Yr1', Elec_bill_withoutDER[0])
    # print('Solar Power : ', sum(solpower))
    average_monthlycashflow = average_annualcashflow / 12
    # print('average_monthlycashflow:', average_monthlycashflow)

    # Calculation of payback period
    payback_year = 0
    for i in range (0,26):
        if (cCF_t[i] < 0 and cCF_t[i+1]>0):
            payback_year = ((i) + (-cCF_t[i]/CF[i+1]))+1
    # print('The payback year is:', payback_year)

    # payback_months = eq_amount / average_monthlycashflow
    # # print('paybackmonths:', payback_months)
    # payback_year = payback_months / 12
    # # print('payback_year:', payback_year)
    # payback_year = payback_year
    # # print('paybackyear:', payback_year)
    total_savings_bill = sum(total_savings)
    roi = ((cum_cashflow) / sum(total_cost)) * 100 * (1 / 25)
    roi = roi
    # roi = ((sum(total_savings)) /sum(total_cost)) * 100 * (1/25)
    # print('Total cost = ' + str(sum(total_cost)))
    # roi = (((1 + (roi_1)) ** (1/25)) - 1) * 100
    return npv, payback_year, cum_cashflow, roi, total_savings_bill, bau_npv, dis_saving, NPV_to_Savings, amount_invested

# print('The NPV is :',financial_calc(x1))

#Define Input Method
inputmethod='customize'

# Optimisation using Scipy

if inputmethod=='optimize':
# SLSQP optimization
    # define objective: maximize Return on investment
    def objective(x1, sign=-1):
        # print(type(x1))
        return sign * financial_calc(x1)[7]


# SHGO optimization
    # define eggholder: maximize saving percentage
    # def eggholder(x1, sign=-1):
    #     # print(type(x1))
    #     return sign * financial_calc(x1)[7]


    # def constraint2(x1):
    #     return 25 - financial_calc(x1)[1]
    #
    # def constraint3(x1):
    #     return financial_calc(x1)[7] - 10
    #
    #
    # con2 = {'type': 'ineq', 'fun': constraint2}
    # con3 = {'type': 'eq', 'fun': constraint3}
    #
    # cons = ([con3])

#Define bounds
    xx = np.zeros(2, dtype=float)
    # print(type(xx))
    if solar & battery:
        xx[0] = 1
        xx[1] = 1
        b = (1, sload)
        b1 = (0, sload * 2)
        bnds = (b, b1)
    else:
        xx[0] = 1
        xx[1] = 0
        b = (1, sload)
        b1 = (0, 0)
        bnds = (b, b1)

# show initial objective
    # SLSQP optimization
    print('Initial Objective: ' + str(-objective(xx)))


    # SHGO optimization
    # print('Initial Objective: ' + str(-eggholder(xx)))

    # Basin Hopping Optimisation
    # print('Initial Objective: ' + str(-objective(xx)))

    # Record value of objective
    # results_shgo = []
    # results_shgo.append(objective(xx))
    # x_list = []
    # x_list.append(xx[0])
    # y_list = []
    # y_list.append(xx[1])
 # For Call back
    def CB(xx):
        # print('In call back')

        x_value = xx[0]
        x_new = x_value
        y_value = xx[1]
        y_new = y_value
        r_value = -objective(xx)
        r_new = r_value
        group1 = (x_new, y_new, r_new)
        # df1 = pd.DataFrame()
        # df1['Solar'] = x_new
        # df1['Battery'] = y_new
        # df1['Func'] = r_new
        # print(df1)
        print(group1)

        return group1
        # return x_new, y_new, r_new,  # df1

# Optimisation function
    # SLSQP optimization
    solution = minimize(objective, xx, method='SLSQP', bounds=bnds, callback=CB, #constraints=cons,
                        options={'ftol': 1e-2, 'disp': True, 'maxiter': 15})

    # SHGO optimization
    # solution1 = shgo(eggholder, bounds=bnds, n=2, iters=1, callback=CB,sampling_method='simplicial')

    # Basin Hopping Optimisation
    # solution2 = basinhopping(objective, xx, niter=10, T=2.0, stepsize=10,
    #                          minimizer_kwargs={'method':'L-BFGS-B'}, take_step=None, accept_test=None,
    #                          callback=CB(xx), interval=20, disp=True, niter_success=None,
    #                          seed=None)

    # Goal seek optiimization
    # roi = 10
    # goal = roi
    #
    # solution = GoalSeek(objective,goal,xx)

# Saving Optimised result
    # SLSQP optimization
    x = solution.x

    # SHGO optimization
    # x = solution1.x
    # xmore = solution1.xl

    # Basin Hopping optimization
    # x = solution2.x

# show final objective
    # SLSQP optimization
    print(solution.x)

    # SHGO optimization
    # print(solution1.x)
    # print(solution1.x1)

    # Basin Hopping optimization
    # print(solution2.x)

# Display Final objective
    # SLSQP optimization
    # print('Final Objective: ' + str(-objective(x)))

    # SHGO optimization
    # print('Final Objective: ' + str(-eggholder(x)))

    # Basin Hopping optimization
    print('Final Objective: ' + str(-objective(x)))

# Printing Solution
    # print solution
    print('Solution')
    print('x1 = ' + str(x[0]))
    print('x2 = ' + str(x[1]))
    # print('x3 = ' + str(x[2]))


    npv1, pback_year, cum_cashflow, roi1, tot_savings, bau_npv, dis_saving, NPV_to_Savings, amount_invested = financial_calc(
        solution.x)

    # end time
    end = time.time()

    # total time taken
    print(f"Runtime of the program is {end - start}")
    runtime = (end - start)
    print('npv = ' + str(npv1))
    print('payback year = ' + str(pback_year))
    print('cumulative cash flow = ' + str(cum_cashflow))
    print('return on investment = ' + str(roi1))
    print('total savings on bill = ' + str(tot_savings))
    print('NPV for BAU = ' + str(bau_npv))
    print('Discounted Total Savings = ' + str(dis_saving))
    print('NPV(BAU) to dis.Savings = ' + str(NPV_to_Savings))
    a1 = [x[0], x[1], npv1, pback_year, cum_cashflow, roi1, tot_savings, runtime, amount_invested]
    print(a1)
    df = pd.DataFrame(a1)
    df = df.round(decimals=3)
    df.index = ['Solar capacity (kW)', 'Storage capacity (kWh)', 'npv of cashflow (INR)', 'payback year(years)',
                'cumulative cash flow(INR)', 'return on investment(%)', 'total savings on bill', 'runtime (sec)',
                'Installation cost (INR)']

    # df.to_excel('scipy_optresult_maxroi_scenario3.xlsx', header=None)

else:
# Customise run inputs and results display
    a=np.zeros(2, dtype=float)
    a[0] = 5 # user input solar capacity
    # a[0]= 26.126338006179388 # user input solar capacity
    a[1]= 0 # user input storage capacity
    npv, payback_year, cum_cashflow, roi, total_savings_bill, bau_npv, dis_saving, NPV_to_Savings, amount_invested=financial_calc(a)

    end = time.time()
    # total time taken
    print(f"Runtime of the program is {end - start}")
    runtime = (end - start)


    print('NPV = ' + str(npv))
    print('payback year = ' + str(payback_year))
    print('cumulative cash flow = ' + str(cum_cashflow))
    print('return on investment = ' + str(roi))
    print('total savings on bill = ' + str(total_savings_bill))
    print('NPV for BAU = ' + str(bau_npv))
    print('Discounted Total Savings = ' + str(dis_saving))
    print('NPV(BAU) to dis.Savings = ' + str(NPV_to_Savings))
    a1 = [a[0], a[1], npv, payback_year, cum_cashflow, roi, total_savings_bill, bau_npv, dis_saving, NPV_to_Savings, amount_invested]
    print(a1)

# end time
end = time.time()

runtime = (end - start)
print('The runtime is:', runtime)
# Load import from csv file
# load = pd.read_csv("load5000(20,20,40,20).csv", header=None)  # need to be replaced with hourly user consumption data for an year
# print('Load:', load)

