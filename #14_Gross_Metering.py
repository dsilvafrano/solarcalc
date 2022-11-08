# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.1
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Package, Input, SQL conn, initial calc., load (& cost) escalation & Elec w/o bill

# %%
# Package, Input & SQL conn
import numpy as np
import pandas as pd
from scipy.optimize import minimize, basinhopping
import sqlite3
import logger
import numpy_financial as npf
import time
## from shgo import shgo 
import requests
from datetime import date, datetime, timedelta
from calendar import monthrange

## from WhatIfAnalysis import GoalSeek
import json

#Inputs from Solsavi
sload = 120 # from user
tariff_id = 2 # consumer type
voltage_id = 2 # voltage type
voltage = "HT"
tariff = 'Industrial'
residence_type = 'Independent House'
metering_id = 8
state_id = 1
state = "Tamil Nadu"
metering_type = "Gross Metering" #Gross Metering, Net Feed In, Net Metering
weekend_consumption_change = 0.5
weekend_consumption_separate = 1 # if this is 1, it means there is weekend consumption
load_input_type = "average_monthly"# month_wise & average_monthly
avg_monthly = 9013.5

# in case of monthwise load_input_type
mc1 = 10232
mc2 = 10252
mc3 = 10124
mc4 = 12568
mc5 = 8164
mc6 = 8108
mc7 = 11944
mc8 = 12636
mc9 = 13560
mc10 = 13396
mc11 = 9812
mc12 = 8880

# Define load distribution for 24hrs weekday
weekday_consumption_6to10 = 10
weekday_consumption_10to18 = 70
weekday_consumption_18to22 = 10
weekday_consumption_22to6 = 10


# Define load distribution for 24hrs weekend


weekend_consumption_6to10 = 20
weekend_consumption_10to18 = 50
weekend_consumption_18to22 = 20
weekend_consumption_22to6 = 10

tou_select = 2
nyr = 26
solar = True
battery = True

x1 = np.zeros(2, dtype=float)
x1[0] = 1 # user input solar capacity
x1[1] = 0 # user input storage capacity

der_deg = 0.01  # solar degradation
bat_type = 1  # battery type = 1 fo Li ion, 0 for lead acid
socin = 0.6
socmax = 0.9
#Ma
batstatus = pd.read_csv("dispatch_strategy(Teddy).csv", header=None)  # need to be replaced with the actual dispatch strategy
# print(batstatus[24:48])

solarpv_subsidy = 0

# starting time
start = time.time()

# fetch solar data using latitude and longitude (get the latitude and longitude using the pincode data
latitude = 9.894568695410525
longitude = 78.07983441534354
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

# Connect to DB and queries from DB
#fetch the dc ac ratio from database

dc_ac_ratio_q = pd.read_sql_query(
    "select value from assumptions_pvwatts where parameter = 'dc_ac_ratio'", conn)
pysam_dc_ac_ratio = float(dc_ac_ratio_q.values[0])

#fetch the inverter replacement year from database

inv_replace_year_q = pd.read_sql_query(
    "select value from assumptions_pvwatts where parameter = 'inverter_replacement_year'", conn)
inv_replace_year = int(inv_replace_year_q.values[0])

#find the battery replacement year based on the lifecycle of batteries

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
# print('DOD :', dod)

socmin = 1 - dod
# print('SOC min', socmin)
battery_replace_year = int(battery_cycle_v / 365)
numb_battery_replacement = int(nyr / (battery_replace_year + 1))
rep_yrs_battery = [i * (battery_replace_year + 1) for i in range(1, numb_battery_replacement + 1)]
rep_yrs_inverter = [inv_replace_year]

#loan related data
##debt fracton fetching from database

debt_fraction_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'debt_fraction'", conn)
# logger.debug("pysam_debt_fraction ="+str((debt_fraction_q.values[0])))
pysam_debt_fraction = float(debt_fraction_q.values[0])

#Loan Rate from database

loan_rate_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'loan_rate'", conn)
# logger.debug(random_no+" "+"loan_rate_q ="+str((loan_rate_q.values[0])))
loan_rate = float(loan_rate_q.values[0]) / 100

#loan period

loan_period_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'loan_term'", conn)
# logger.debug(random_no+" "+"loan_rate_q ="+str((loan_rate_q.values[0])))
loan_period = float(loan_period_q.values[0])

#cost escalation/inflation rate

cost_esc_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'inflation_rate'", conn)
# logger.debug(random_no+" "+"cost_esc_q ="+str((lcost_esc_q.values[0])))
cost_esc = cost_esc_q.values[0] / 100
# print('Cost escalation :', cost_esc)

#get discount rate for calculating npv

dis_factor_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'real_discount_rate_show'",conn)
dis_factor = dis_factor_q.values[0] / 100
# print('Discount rate :', dis_factor)

#cost escalation/inflation rate

cost_esc_q = pd.read_sql_query("select value from assumptions_cashloan where parameter = 'inflation_rate'", conn)
cost_esc = cost_esc_q.values[0] / 100
# print('Cost escalation :', cost_esc)

#load escalation

load_esc_q = pd.read_sql_query("select value from assumptions_grid where parameter = 'load_escalation'", conn)
load_esc = load_esc_q.values[0] / 100
# print('Load escalation:', load_esc)


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

#Retrieve financial suggestions

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

#     print('The capital cost is:', capital_cost)

    return max_limit_pv, min_limit_pv, capital_cost, inverter_cost, hybrid_inverter, pv_cost, li_ion, lead_acid

# Replacement cost calculation

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
print('The replacement costs are:',sum(replacement_cost(10, 0)[0]), sum(replacement_cost(10, 0)[1]) )
    
# Calculate the hourly load from avg_monthly or monthwise input

## Percentage of consumption in weekday & weekend
user_load = np.array(([0] * 8760), dtype=float)

n = len(user_load)  # number of load values



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
    
# Populating the days of the month
## Populating the days of the month
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

        # print("user load full matrix ",user_load)
        # print(len(user_load))
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

    print(sum(Year1_monthly))
    # print(len(Year1_monthly))

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

            days_in_month = monthrange(2021, i)[1]
            print(days_in_month)
            given_date = datetime(year=2021, month=i, day=1).date()
#             print(given_date)
            first_day_of_month = given_date.replace(day=1)
            weekdays_in_month = np.busday_count(first_day_of_month, last_day_of_month(given_date))
#             print(last_day_of_month(given_date))
            weekends_in_month = days_in_month - weekdays_in_month

            print("weekdays ", weekdays_in_month)
            print("weekends ", weekends_in_month)

            consumption_weekday = mc[i - 1] / (
                        weekdays_in_month + (weekends_in_month * (weekend_consumption_change + 1)))
            consumption_weekend = ((weekend_consumption_change + 1) * consumption_weekday)
            print('consumption_weekday ', consumption_weekday)

            print('consumption_weekend ', consumption_weekend)

#             print(weekend_consumption_6to10n)
#             print(weekend_consumption_10to18n)
#             print(weekend_consumption_18to22n)
#             print(weekend_consumption_22to6n)

            weekday_24 = np.array(([0] * 24), dtype=float)
            weekend_24 = np.array(([0] * 24), dtype=float)

            weekday_24[6:10] = round((weekday_consumption_6to10n * consumption_weekday / 4), 3)
            weekday_24[10:18] = round((weekday_consumption_10to18n * consumption_weekday / 8), 3)
            weekday_24[18:22] = round((weekday_consumption_18to22n * consumption_weekday / 4), 3)
            weekday_24[22:24] = round((weekday_consumption_22to6n * consumption_weekday * 0.25 / 2), 3)
            weekday_24[0:6] = round((weekday_consumption_22to6n * consumption_weekday * 0.75 / 6), 3)
            print(weekday_24)
            
            weekend_24[6:10] = round((weekend_consumption_6to10n * consumption_weekend / 4), 3)
            weekend_24[10:18] = round((weekend_consumption_10to18n * consumption_weekend / 8), 3)
            weekend_24[18:22] = round((weekend_consumption_18to22n * consumption_weekend / 4), 3)
            weekend_24[22:24] = round((weekend_consumption_22to6n * consumption_weekend * 0.25 / 2), 3)
            weekend_24[0:6] = round((weekend_consumption_22to6n * consumption_weekend * 0.75 / 6), 3)
            print(weekend_24)
#             year = date.today().year

            begin_year = date(2021, 1, 1)
            end_year = date(2021, 12, 31)
            one_day = timedelta(days=1)

            # load_value = [];

            # print(begin_year, end_year, one_day)
            # Allocate number of hrs for each month
            # Allocate the beginning of each month
            B = [0,744,1416,2160,2880,3624,4344,5088,5832,6552,7296,8016,8760]
            E = [744,672,744,720,744,720,744,744,720,744,720,744]

            a = 0
            b = 24
            user_load_t = np.array(([0] * E[i-1]), dtype=float)
            next_day = first_day_of_month

            
            for day in range(0, days_in_month):  # includes potential leap year
                if next_day > last_day_of_month(given_date):
                    break
                if (next_day.weekday() > 4):
                    user_load_t[a:b] = weekend_24[0:24]
                    a = b
                    b = b + 24
                else:
                    user_load_t[a:b] = weekday_24[0:24]
                    a = b
                    b = b + 24
                next_day += one_day
                user_load[B[i-1]:B[i]] = user_load_t[0:E[i-1]]
    print('user_load', sum(user_load))
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

    print((Year1_monthly))
    # print(len(Year1_monthly))
    # print('Load :', sum(user_load))
    # print(len(user_load))
    
# Apply load escalation annually
def annual_load_escalation(load_esc):
    annual_load = []
    #Applying escalation for 25 years
    for n in range(0,26):
        esc_load_n = user_load
        esc_load_n = esc_load_n * (1 + (n * load_esc))
        annual_load.append(esc_load_n)

    return annual_load

# print('The annual load for year:', sum(annual_load_escalation(load_esc)[0]))

# Find monthly load from escalated load
## Each month cumulative has been calculated with applied escalation
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

# print('The monthly escalated load:',(monthly_cum_load(0))[8759])

# Energy and Fixed charge calculation
## Retrieve the max applicable slab and tiers in EC
slab_id_q = pd.read_sql_query("select slab_id, period, tier, min, maximum, energy_charge from slabs_mapping where state_id =" 
                              + str(state_id) + " and tarriff_type_id = " + str(tariff_id) + " and voltage_type_id = " 
                              + str(voltage_id) + " and metering_type_id = " + str(metering_id),conn)
slab_id_t = slab_id_q
# print('The EC table is:', slab_id_t)
slab_id_m = slab_id_t['slab_id'].max()
# print('The max slab is', slab_id_m)
tier_m = slab_id_t['tier'].max()
# print('The tier slab is', tier_m)

# Retrieve and check the fixed charge calculation type
charge_calculation_q = pd.read_sql_query("select charge_calculation from fixedcharge where tariff_id=" + str(tariff_id) 
                                         + " and voltage_id=" + str(voltage_id) + " and state_id=" + str(state_id), conn)
charge_calculation = int(charge_calculation_q.values[0])
if int(charge_calculation) == 0:
    fixed_charge = 0
elif charge_calculation == 1:
    fixed_charges_q = pd.read_sql_query("select fixed_charge from fixedcharge where state_id =" + str(state_id) 
                                        + " and voltage_id = " + str(voltage_id),conn)
    fixed_charge = float(fixed_charges_q.values[0])
elif charge_calculation == 2:
    
# Check if the fixed charge table has tiers
    tier_check_q = pd.read_sql_query("select tier from fixedcharge where tariff_id=" + str(tariff_id) + " and voltage_id=" 
                                     + str(voltage_id) + " and state_id=" + str(state_id), conn)
    tier_check = tier_check_q
    # print('Tier check:', tier_check)
    max_tier_check = tier_check['tier'].max()
    # print('Max tier check:', max_tier_check)
    min_tier_check = tier_check['tier'].min()
    # print('Min tier check:', min_tier_check)
    
# Retrieve fixed charge from DB and calculation
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

# Fixed charge and Energy charge escalation
# Fixed charge escalation
# apply escalation to the fixed charge calculation
def fixed_charge_esc(n):
    FC = fixed_charge_h

    for j in range(n, n+1):
        FC_esc = FC * (1 + (j * cost_esc))

    # print('The escalated fixed charge is:', FC_esc)

    return FC_esc
# print('The fixed charge is:', fixed_charge_esc(0))

# Energy charge escalation (normal)
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

# Energy charge escalation (peak)
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

#Energy charge escalation (off peak)
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

# Slab selection for bill without system
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

# Applicable TOU charges for the user
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

# Selection of Energy charge
## Build a cost matrix with cost for each hour (8760)
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

# %% [markdown]
# # Electricity bill without system

# %%
# Calculate the bill without system

def Elec_bill_wo_system(n):
    # starting time
    # start1 = time.time()
    Elec_bill_wo_system = 0
    bill_load = annual_load_escalation(load_esc)[n]
    bill_cost = charge_selection_n(n)
    fixed_charge_h_n = fixed_charge_esc(n)
    bill_annual = 0

    for i in range(0, 8760):
        bill_annual = bill_annual + (fixed_charge_h_n + (bill_load[i] * bill_cost[i]))

    Elec_bill_wo_system = bill_annual

    # print('The annual bill w/o system is:', Elec_bill_wo_system)
    # end time
    # end1 = time.time()
    # runtime1 = (end1 - start1)
    # print('The runtime of elec bill wo system is:', runtime1)

    return Elec_bill_wo_system

print('The annual bill is:', Elec_bill_wo_system(0))


# %% [markdown]
# # Power balance of the system

# %%
# power balance for the selected size
def power_balance(x1, g):
    # starting time
    start2 = time.time()
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

    # end time
    end2 = time.time()
    runtime2 = (end2 - start2)
    # print('The runtime for power balance is:',runtime2)
#     print('Total solar :', sum_solar)
    # print('Export:', sum(excessder))
    return battery_power, gridpower, excessder, solar_power, load, df1, sum_battery, sum_grid, sum_export, sum_solar, sum_solar_load


print('The power balance details are:', (power_balance(x1, 0)[5]))

print('The total solar generation is:', power_balance(x1, 0)[9])
print('The total solar contribution is:', power_balance(x1, 0)[10])
print('The total battery contribution is:', power_balance(x1, 0)[6])
print('The total grid contribution is:', power_balance(x1, 0)[7])
print('The total export contribution is:', power_balance(x1, 0)[8])


# %% [markdown]
# # Apply escalation to grid and export, also save in monthly for 25 years

# %%
def grid_monthly_cum_load(x1, n):
    # starting time
    start3 = time.time()
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
    # end time
    end3 = time.time()
    runtime3 = (end3 - start3)
    print('The runtime for grid monthly cum. load is:', runtime3)

    return g_monthly_cum_load, e_monthly_cum_load
# print('The cumulative export is:', (grid_monthly_cum_load(0)[1]))

# %% [markdown]
# # Slab selection after system, based on avg monthly consumption

# %%
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


# %% [markdown]
# # Apply escalation to the energy charge for 25 years

# %%
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

# %% [markdown]
# # Selection of applicable charges

# %% [markdown]
# ## Build a cost matrix with 8760 points

# %%
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

# %% [markdown]
# ## Network charges and escalation

# %%
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


# %% [markdown]
# # Electricity bill with System

# %%
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
            nf_bill_annual = nf_bill_annual + (fixed_charge_h_n + (((cons_grid[g] * cost_grid[g]) - (cons_export[g] * cost_solar)) + (cons_solar[g] * network_charge)))

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

print(' The bill with system is:', Elec_bill_w_sys(x1, 0))
