import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt

cwd = os.getcwd()
datafile = os.path.join(cwd, "data/respiratory_data.csv")

df = pd.read_csv(datafile, thousands=',', usecols=["Week Ending Date", "Total COVID-19 Admissions", "Total Influenza Admissions", "Geographic aggregation"])
df["Week Ending Date"] = pd.to_datetime(df["Week Ending Date"])
omicron_time = df[df["Week Ending Date"] <= pd.to_datetime("03/01/2022")]
omicron_time = omicron_time[omicron_time["Geographic aggregation"]=="USA"]

plt.plot(range(0,82), omicron_time["Total COVID-19 Admissions"])
plt.plot(range(0,82), omicron_time["Total Influenza Admissions"])
plt.show()