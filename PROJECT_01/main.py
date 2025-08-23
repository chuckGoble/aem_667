"""AEM 667: Project 1 Code"""
#Imports
from glob import glob
from matplotlib import pyplot as plt
from pathlib import Path
import os
import sys
import numpy as np
import pandas as pd

# SETUP Path Constants
CODE = Path(os.getcwd())

# WGS-84 Constants
MU = 3.986005 * 10**14  # m^3/s^2
SPEED_OF_LIGHT = 2.99792458 * 10**8  # m/s
EARTH_ROT_RATE = 7.2921151467 * 10**-5  # rad/s
PI = 3.1415926535898
FLATTENING = 1 / 298.257223563
ECCENTRICITY = np.sqrt(FLATTENING * (2 - FLATTENING))
EQ_RAD = 6378137  # m
POL_RAD = EQ_RAD * (1 - FLATTENING)


def make_icp_df(path_list):
    sats = list()
    for file_path in path_list:
        df = pd.read_table(file_path, names=list(range(8)))
        df = df.loc[:, [0, 1, 5, 6, 7]]
        df.columns = ['TIME_SEC', 'CARRIER_PHASE_CYCLES', 'CYCLE_SLIP_COUNTER',
                     'SNR', 'PSEUDORANGE_M']
        df['SAT'] = Path(file_path).name.split('_')[1].split('.')[0].split('t')[1]
        sats.append(df)
    return pd.concat(sats)


if __name__ == '__main__':
    # Load Base Station ECEF Data
    base_ecef = glob(str(CODE / '*' / '*ecef*.txt'))
    df_base_ecef = pd.read_table(base_ecef[0], names=list(range(18)))
    df_base_ecef = df_base_ecef.loc[:, list(range(8))].copy()
    df_base_ecef.columns = ['TIME_SEC', 'GPS_WEEK', 
                 'ECEF_POS_X_M', 'ECEF_POS_Y_M', 'ECEF_POS_Z_M',
                 'ECEF_VEL_X_MPS', 'ECEF_VEL_Y_MPS', 'ECEF_VEL_Z_MPS']

    # Load Base Station Rx Data 
    base_icp = glob(str(CODE / '*base*' / '*icp*.txt'))
    df_base = make_icp_df(base_icp)
    base_unique_sats = df_base['SAT'].unique().tolist()

    # Load Rover Rx Data
    rvr_icp = glob(str(CODE / '*rover*' / '*icp*.txt'))
    df_rvr = make_icp_df(rvr_icp)
    rvr_unique_sats = df_rvr['SAT'].unique().tolist()

    # Load Broadcast Data
    sats = set(base_unique_sats + rvr_unique_sats)
    gps_eph = glob(str(CODE / '*' / '*broadcast*.csv'))
    df_eph = pd.read_csv(gps_eph[0])
    gps_sat_state = glob(str(CODE / '*' / '*sat_state*.csv'))
    df_gps = pd.read_csv(gps_sat_state[0])

    # Identify all the unique Satellites between the Base and Rover
    unique_sats = list(set(base_unique_sats + rvr_unique_sats))

    # Task 001: Compute average Base receiver position
    base_stats_full_data = df_base_ecef.describe()
    df_base_ecef = df_base_ecef.loc[18:, :].reset_index(drop=True)
    base_stats_full_data = df_base_ecef.describe()
    base_stats = df_base_ecef.describe()
    base_origin_ecef = np.array([[base_stats['ECEF_POS_X_M']['mean']],
                                  [base_stats['ECEF_POS_Y_M']['mean']],
                                   [base_stats['ECEF_POS_Z_M']['mean']]])
    
    plt.figure()
    plt.subplot(3, 1, 1)
    plt.scatter(df_base_ecef['TIME_SEC'], df_base_ecef['ECEF_POS_X_M'])
    plt.plot(df_base_ecef['TIME_SEC'],
             base_stats['ECEF_POS_X_M']['mean'] * np.ones(len(df_base_ecef['TIME_SEC'])), color='red')
    plt.grid()
    plt.title('GPS WEEK:' + ' ' + str(df_base_ecef['GPS_WEEK'][0]))
    plt.xlabel('GPS Time of Week [sec]')
    plt.ylabel('ECEF X-Axis [m]')
    plt.legend(['Data', 'Mean'])

    plt.subplot(3, 1, 2)
    plt.scatter(df_base_ecef['TIME_SEC'], df_base_ecef['ECEF_POS_Y_M'])
    plt.plot(df_base_ecef['TIME_SEC'],
             base_stats['ECEF_POS_Y_M']['mean'] * np.ones(len(df_base_ecef['TIME_SEC'])), color='red')
    plt.grid()
    plt.xlabel('GPS Time of Week [sec]')
    plt.ylabel('ECEF Y-Axis [m]')
    plt.legend(['Data', 'Mean'])

    plt.subplot(3, 1, 3)
    plt.scatter(df_base_ecef['TIME_SEC'], df_base_ecef['ECEF_POS_Z_M'])
    plt.plot(df_base_ecef['TIME_SEC'],
             base_stats['ECEF_POS_Z_M']['mean'] * np.ones(len(df_base_ecef['TIME_SEC'])), color='red')
    plt.grid()
    plt.xlabel('GPS Time of Week [sec]')
    plt.ylabel('ECEF Z-Axis [m]')
    plt.legend(['Data', 'Mean'])

    # Task 002: Compute Line-Of-Sight (LOS) vectors for the Base and Rover in NED coordinates
    # Load the PRN State Data
    prn_state_list = glob(str(CODE / '*' / '*gps*state*prn*.csv'))
    prn_state_list = [pd.read_csv(ii) for ii in prn_state_list]
    prn_states = pd.concat(prn_state_list)

    # Compute the relative position vector between the orbit and base/rover
    base_los = prn_states[:, 'SAT_POS_X_m':'SAT_POS_Z_M'] - base_origin_ecef.T
    rover_los = prn_states[:, 'SAT_POS_X_m':'SAT_POS_Z_M'] - base_origin_ecef.T

    # TODO: Plots