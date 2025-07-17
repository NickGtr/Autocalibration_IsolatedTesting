import pandas as pd
import numpy as np
import matplotlib.pyplot as plt
from matplotlib.ticker import PercentFormatter
import math
import os
from scipy.stats import lognorm

base_path = "/Data/gautt/bundle_adjustment/autocalibration_project/results"

def QRvsIAC_10vs50(num_reconstructions, translation_range, ylim = None):
    # Load both datasets
    file_name_10 = f"QRvsIAC_10Cams_{num_reconstructions}Reconstr_{translation_range}Transl.txt"
    file_name_50 = f"QRvsIAC_50Cams_{num_reconstructions}Reconstr_{translation_range}Transl.txt"
    path_10 = os.path.join(base_path, file_name_10)
    path_50 = os.path.join(base_path, file_name_50)
    df_10 = pd.read_csv(path_10, sep=r'\s+')
    df_50 = pd.read_csv(path_50, sep=r'\s+')

    # Create side-by-side plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    # Plot for 10 cameras
    axes[0].plot(df_10['FocalLength'], df_10['IACloss'], marker='s', label='IACloss')
    axes[0].plot(df_10['FocalLength'], df_10['QRloss'], marker='o', label='QRloss')
    axes[0].set_title('10 Cameras')
    axes[0].set_xlabel('Focal Length (in mm)')
    axes[0].set_ylabel('Losses (in mm)')
    axes[0].legend()
    axes[0].grid(True)

    # Plot for 50 cameras
    axes[1].plot(df_50['FocalLength'], df_50['IACloss'], marker='s', label='IACloss')
    axes[1].plot(df_50['FocalLength'], df_50['QRloss'], marker='o', label='QRloss')
    axes[1].set_title('50 Cameras')
    axes[1].set_xlabel('Focal Length (in mm)')
    axes[1].legend()
    axes[1].grid(True)
    axes[1].tick_params(labelleft=True)  # <--- This line ensures y-axis labels are shown

    # Main title and layout
    fig.suptitle(f'FIXED FOCAL LENGHT - 2004 paper - QR vs IAC Loss depending on Equivalent Focal Length : {num_reconstructions} reconstructions per focal length; {translation_range} translation_range')
    if (ylim is not None):
        plt.ylim(0,ylim)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()

def QRvsIAC_10vs50_1997paper(num_reconstructions, translation_range, ylim = None):
    # Load both datasets
    file_name_10 = f"QRvsIAC_10Cams_{num_reconstructions}Reconstr_{translation_range}Transl_1997paper.txt"
    file_name_50 = f"QRvsIAC_50Cams_{num_reconstructions}Reconstr_{translation_range}Transl_1997paper.txt"
    path_10 = os.path.join(base_path, file_name_10)
    path_50 = os.path.join(base_path, file_name_50)
    df_10 = pd.read_csv(path_10, sep=r'\s+')
    df_50 = pd.read_csv(path_50, sep=r'\s+')

    # Create side-by-side plots
    fig, axes = plt.subplots(1, 2, figsize=(14, 5), sharey=True)

    # Plot for 10 cameras
    axes[0].plot(df_10['FocalLength'], df_10['IACloss'], marker='s', label='IACloss')
    axes[0].plot(df_10['FocalLength'], df_10['QRloss'], marker='o', label='QRloss')
    axes[0].set_title('10 Cameras')
    axes[0].set_xlabel('Focal Length (in mm)')
    axes[0].set_ylabel('Losses (in mm)')
    axes[0].legend()
    axes[0].grid(True)

    # Plot for 50 cameras
    axes[1].plot(df_50['FocalLength'], df_50['IACloss'], marker='s', label='IACloss')
    axes[1].plot(df_50['FocalLength'], df_50['QRloss'], marker='o', label='QRloss')
    axes[1].set_title('50 Cameras')
    axes[1].set_xlabel('Focal Length (in mm)')
    axes[1].legend()
    axes[1].grid(True)
    axes[1].tick_params(labelleft=True)  # <--- This line ensures y-axis labels are shown

    # Main title and layout
    fig.suptitle(f'FIXED FOCAL LENGHT - 1997 PAPER - QR vs IAC Loss depending on Equivalent Focal Length : {num_reconstructions} reconstructions per focal length; {translation_range} translation_range')
    if (ylim is not None):
        plt.ylim(0,ylim)
    plt.tight_layout(rect=[0, 0, 1, 0.95])
    plt.show()