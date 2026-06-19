"""
Script used to add baseflow to quickflow hydrographs (if using RORB) and/or route the flows through the dam
"""
import sys
import os
import pandas as pd
from lib.Routing import Router

# Get the config information
routing_sheet = sys.argv[1]
print('Opening routing spreadsheet: ', routing_sheet)
routing_df = pd.read_excel(routing_sheet, sheet_name=0)
routing_df = routing_df.loc[routing_df['Include'] == 'yes']
print(routing_df)

# Run through the simulations
for index, sim in routing_df.iterrows():
    Router(sim)
