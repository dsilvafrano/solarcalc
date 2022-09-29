# ---
# jupyter:
#   jupytext:
#     formats: ipynb,py:percent
#     text_representation:
#       extension: .py
#       format_name: percent
#       format_version: '1.3'
#       jupytext_version: 1.14.0
#   kernelspec:
#     display_name: Python 3 (ipykernel)
#     language: python
#     name: python3
# ---

# %% [markdown]
# # Replicate Financial Computation
#
# We will replicate the tier computations in this notebook.
#
# ## User Input Values and Options

# %%
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

user_inputs {
    "battery": True,
    "solar": True,
    "tariff": "Domestic",
    "residence_type": 'Independent House',
}


# %% [markdown]
# ## Find Total Installation Cost

# %%
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

