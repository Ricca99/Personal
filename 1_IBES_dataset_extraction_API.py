# -*- coding: utf-8 -*-
"""
Created on Thu Jul 17 10:37:01 2025

@author: meliric
"""

#######################################################################################
#Code to extract Eikon projections data
#######################################################################################

### Import relevant libraries and refinitiv python api
import pandas as pd
from pandas import isna
import os
import eikon as ek
import numpy as np
from glob import glob
from connectors import devo
import pyodbc
con = pyodbc.connect('DSN=DEVO Impala 64bit',autocommit=True)
ek.set_app_key('a8b11a96a951495f95c52893aebeba888985c789')

### Sample information
# Import dataset
extract_tickers = pd.read_excel(
    "P:/ECB business areas/DGMF/SRF/Data management/Market data/IBES expectations/sample_info.xlsx",
    sheet_name="sample_info",
    usecols="A:O",
    dtype=str)
# get the ric tickers from excel so as to proceed with extraction
tickers = extract_tickers['ric'].dropna().unique().tolist()

### the tickers in the current sample are 88 - cross-checked with the following vector 
# check_tickers = ['BBVA.MC','BAMI.MI','SAN.MC','BNPP.PA','CRDI.MI','MDBI.MI','KBC.BR',
#             'INGA.AS','ERST.VI','ISP.MI','EMII.MI','ABNd.AS','SOGN.PA','BKT.MC',
#             'SABE.MC','CABK.MC','CBKG.DE','DBKGn.DE','BIRG.I','CAGR.PA',
#             'NDASE.ST','BAWG.VI','FBK.MI','RBIV.VI','ARLG.DE','PBBG.DE',
#             'UNI.MC','EURBr.AT','NBGr.AT','BPSI.MI','BCP.LS','EMBI.MI','DANSKE.CO','JYSK.CO','SYDB.CO',
#             'SEBa.ST','SHBa.ST','SWEDa.ST','HSBA.L','BARC.L','NWG.L','LLOY.L','STAN.L','JPM.N',
#             'PNC.N','BAC.N','C.N','CMA.N','CFG.N','RF.N','MTB.N',
#             'FITB.OQ','HBAN.OQ','PBCT.OQ^D22','KEY.N','USB.N',
#             'TFC.N','SBNY.OQ','SIVB.OQ','FRC.N','ZION.OQ','WFC.N','8306.T',
#             '8411.T','8316.T','8308.T','8604.T','8304.T','8331.T','8354.T',
#             '7186.T','8303.T','8355.T','UBSG.S','CSGN.S','ARLG.DE','SAB1L.VL','ADKO.VI',
#             'AIBG.I','8589934341','ACBr.AT','BMED.MI','BMPS.MI',
#             'BOCH.CY','HBNK.CY','NDAFI.HE','NLBR.LJ','BOPr.AT','IL0A.I']

### Extract the data - try to run all together first. IF not working, run separately each chunck
years = ['FY1','FY2','FY3']

# Net income (after tax, in billions)
for j in years: 
    name_ = f"df_NetProfit_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.NetProfitMean.calcdate','TR.NetProfitMean.fperiod','TR.NetProfitMean(Scale=9)',
                                    'TR.NetProfitMedian(Scale=9)','TR.NetProfitEstStdDev(Scale=9)','TR.NetProfitNumOfEstimates',
                                    'TR.NetProfitHigh(Scale=9)','TR.NetProfitLow(Scale=9)'],
                          parameters = {'SDate' : '2025-01-01', 'EDate': '1D', 'Period':j, 'Frq':'C', 'Curn':'EUR'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# ROE (in %)
for j in years: 
    name_ = f"df_ROE_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.ROEMean.calcdate','TR.ROEMean.fperiod','TR.ROEMean','TR.ROEMedian',
                                    'TR.ROEEstStdDev','TR.ROENumOfEstimates','TR.ROEHigh','TR.ROELow'],
                          parameters = {'SDate' : '2025-01-01', 'EDate': '1D', 'Period':j, 'Frq':'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Fees and commission income (in billions)
for j in years:  
    name_ = f"df_FeesandCommissionsInc_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.FeesandCommissionsIncMean.calcdate','TR.FeesandCommissionsIncMean.fperiod',
                                    'TR.FeesandCommissionsIncMean(Scale=9)','TR.FeesandCommissionsIncMedian(Scale=9)','TR.FeesandCommissionsIncStdDev(Scale=9)',
                                    'TR.FeesandCommissionsIncNumOfEstimates','TR.FeesandCommissionsIncHigh(Scale=9)','TR.FeesandCommissionsIncLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

 # Total loans (in billions)
for j in years: 
    name_ = f"df_Loans_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.LoansMean.calcdate','TR.LoansMean.fperiod','TR.LoansMean(Scale=9)',
                                    'TR.LoansMedian(Scale=9)','TR.LoansStdDev(Scale=9)','TR.LoansNumOfEstimates',
                                    'TR.LoansHigh(Scale=9)','TR.LoansLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C', 'Curn':'EUR'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Loss provisions (in billions)
# A measure of the expense to account for future losses on customer loan defaults
for j in years:  
    name_ = f"df_LoanLossProvisions_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.LoanLossProvisionsMean.calcdate','TR.LoanLossProvisionsMean.fperiod','TR.LoanLossProvisionsMean(Scale=9)',
                                    'TR.LoanLossProvisionsMedian(Scale=9)','TR.LoanLossProvisionsStdDev(Scale=9)','TR.LoanLossProvisionsNumOfEstimates',
                                    'TR.LoanLossProvisionsHigh(Scale=9)','TR.LoanLossProvisionsLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Net interest income (in billions)
for j in years:  
    name_ = f"df_NetIntInc_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.NetIntIncMean.calcdate','TR.NetIntIncMean.fperiod','TR.NetIntIncMean(Scale=9)','TR.NetIntIncMedian(Scale=9)',
                                    'TR.NetIntIncStdDev(Scale=9)','TR.NetIntIncNumOfEstimates','TR.NetIntIncHigh(Scale=9)','TR.NetIntIncLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

 # Net trading income (in billions)
for j in years: 
    name_ = f"df_TradingInc_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.TradingIncMean.calcdate','TR.TradingIncMean.fperiod','TR.TradingIncMean(Scale=9)',
                                    'TR.TradingIncMedian(Scale=9)','TR.TradingIncStdDev(Scale=9)','TR.TradingIncNumOfEstimates',
                                    'TR.TradingIncHigh(Scale=9)','TR.TradingIncLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Operating expenses (in billions)
for j in years:  
    name_ = f"df_OperatingExp_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.OperatingExpMean.calcdate','TR.OperatingExpMean.fperiod','TR.OperatingExpMean(Scale=9)',
                                    'TR.OperatingExpMedian(Scale=9)','TR.OperatingExpStdDev(Scale=9)','TR.OperatingExpNumOfEstimates',
                                    'TR.OperatingExpHigh(Scale=9)','TR.OperatingExpLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Operating income (in billions)
for j in years:  
    name_ = f"df_OPR_{j}"  
    df, err = ek.get_data(tickers,
                         fields = ['TR.OPRMean.calcdate','TR.OPRMean.fperiod','TR.OPRMean(Scale=9)','TR.OPRMedian(Scale=9)',
                                   'TR.OPRStdDev(Scale=9)','TR.OPRNumOfEstimates','TR.OPRHigh(Scale=9)','TR.OPRLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

 # Pre-tax profit (in billions)
 # Company's net income before tax expense
for j in years: 
    name_ = f"df_PreTaxProfit_{j}"  
    df, err = ek.get_data(tickers,
                         fields = ['TR.PreTaxProfitMean.calcdate','TR.PreTaxProfitMean.fperiod','TR.PreTaxProfitMean(Scale=9)',
                                   'TR.PreTaxProfitMedian(Scale=9)','TR.PreTaxProfitEstStdDev(Scale=9)','TR.PreTaxProfitNumOfEstimates',
                                   'TR.PreTaxProfitHigh(Scale=9)','TR.PreTaxProfitLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Total Assets (in billions)
for j in years:  
    name_ = f"df_TotalAssets_{j}"
    df, err = ek.get_data(tickers,
                          fields = ['TR.TotalAssetsMean.calcdate','TR.TotalAssetsMean.fperiod','TR.TotalAssetsMean(Scale=9)','TR.TotalAssetsMedian(Scale=9)',
                                    'TR.TotalAssetsStdDev(Scale=9)','TR.TotalAssetsNumOfEstimates','TR.TotalAssetsHigh(Scale=9)','TR.TotalAssetsLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Total Equity (in billions)
# Computed as firm's total assets minus its total liabilities
for j in years:  
    name_ = f"df_ShareholdersEquity_{j}" 
    df, err = ek.get_data(tickers,
                          fields = ['TR.ShareholdersEquityMean.calcdate','TR.ShareholdersEquityMean.fperiod','TR.ShareholdersEquityMean(Scale=9)',
                                    'TR.ShareholdersEquityMedian(Scale=9)','TR.ShareholdersEquityStdDev(Scale=9)','TR.ShareholdersEquityNumOfEstimates',
                                    'TR.ShareholdersEquityHigh(Scale=9)','TR.ShareholdersEquityLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Total income (in billions)
# it is the interest earned plus other incomes before deducting interest expense
for j in years:  
    name_ = f"df_TotalInc_{j}"  
    df, err = ek.get_data(tickers,
                         fields = ['TR.TotalIncMean.calcdate','TR.TotalIncMean.fperiod','TR.TotalIncMean(Scale=9)','TR.TotalIncMedian(Scale=9)',
                                   'TR.TotalIncStdDev(Scale=9)','TR.TotalIncNumOfEstimates','TR.TotalIncHigh(Scale=9)','TR.TotalIncLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Non-interest income (in income)
for j in years:  
    name_ = f"df_TotalnonIntRev_{j}"  
    df, err = ek.get_data(tickers, 
                          fields = ['TR.TotalnonIntRevMean.calcdate','TR.TotalnonIntRevMean.fperiod','TR.TotalnonIntRevMean(Scale=9)',
                                    'TR.TotalnonIntRevMedian(Scale=9)','TR.TotalnonIntRevStdDev(Scale=9)','TR.TotalnonIntRevNumOfEstimates',
                                   'TR.TotalnonIntRevHigh(Scale=9)','TR.TotalnonIntRevLow(Scale=9)'],
                          parameters={'SDate': '2025-01-01', 'EDate': '1D', 'Period': j, 'Frq': 'C'})
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

 # Earnings per Share
for j in years:
    name_ = f"df_EPS_{j}" 
    df, err = ek.get_data(tickers,
                          fields = ['TR.EPSMean.calcdate','TR.EPSMean.fperiod','TR.EPSMean','TR.EPSMedian','TR.EPSStdDev',
                                    'TR.EPSNumOfEstimates','TR.EPSHigh','TR.EPSLow'],
                          parameters = {'SDate' : '2025-01-01', 'EDate': '1D', 'Period':j, 'Frq':'C'})    
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Dividend per Share
for j in years: 
    name_ = f"df_DPS_{j}" 
    df, err = ek.get_data(tickers,
                          fields = ['TR.DPSMean.calcdate','TR.DPSMean.fperiod','TR.DPSMean','TR.DPSMedian','TR.DPSStdDev',
                                    'TR.DPSNumOfEstimates','TR.DPSHigh','TR.DPSLow'],
                          parameters = {'SDate' : '2025-01-01', 'EDate': '1D', 'Period':j, 'Frq':'C'})    
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

# Dividend per Share Yield
for j in years: 
    name_ = f"df_DPSYield_{j}" 
    df, err = ek.get_data(tickers,
                          fields = ['TR.DPSMeanYield.calcdate','TR.DPSMeanYield.fperiod','TR.DPSMeanYield','TR.DPSMeanYieldMedian',
                                    'TR.DPSMeanYieldStdDev','TR.DPSMeanYieldNumOfEstimates','TR.DPSMeanYieldHigh','TR.DPSMeanYieldLow'],
                          parameters = {'SDate' : '2025-01-01', 'EDate': '1D', 'Period':j, 'Frq':'C'})    
    df = df.rename(columns={'Calc Date': 'Date'})
    globals()[name_] = df
    continue

#######################################################################################
# Export the dataframes into parquete files in local folder - not necessary to run
#######################################################################################

# # Ensure the folder exists
# save_path = r"P:\ECB business areas\DGMF\SRF\Data management\Market data\IBES expectations\new procedure_in progress\temp datasets"
# os.makedirs(save_path, exist_ok=True)

# # Dataset base names
# dataset_names = [
#     "df_NetProfit",
#     "df_ROE",
#     "df_FeesandCommissionsInc",
#     "df_Loans",
#     "df_LoanLossProvisions",
#     "df_NetIntInc",
#     "df_TradingInc",
#     "df_OperatingExp",
#     "df_PreTaxProfit",
#     "df_OPR",
#     "df_TotalAssets",
#     "df_TotalInc",
#     "df_ShareholdersEquity",
#     "df_TotalnonIntRev",
#     "df_EPS",
#     "df_DPS",
#     "df_DPSYield"
# ]

# # Loop and save
# for j in years:
#     for base_name in dataset_names:
#         var_name = f"{base_name}_{j}"
#         try:
#             df = globals()[var_name]   # Get dataframe by name
#             file_path = os.path.join(save_path, f"{var_name}.parquet")
#             df.to_parquet(file_path, index=False)
#             print(f"Exported {var_name} -> {file_path}")
#         except KeyError:
#             print(f" {var_name} not found in globals()")
#         except Exception as e:
#             print(f" Error exporting {var_name}: {e}")

# # Dictionary to store loaded DataFrames
# dataframes = {}

# for j in years:
#     for base_name in dataset_names:
#         var_name = f"{base_name}_{j}"
#         file_path = os.path.join(save_path, f"{var_name}.parquet")

#         if os.path.exists(file_path):
#             df = pd.read_parquet(file_path)
#             dataframes[var_name] = df   # store in dictionary
#             globals()[var_name] = df    # optional: recreate as variable in namespace
#             print(f"Imported {var_name}")
#         else:
#             print(f" File not found: {file_path}")
            
###############################################################################
# DATA CLEANING 
###############################################################################

### (1) check for missing values and save files in the folder 'missing'
# missing_path = r"P:\ECB business areas\DGMF\SRF\Data management\Market data\IBES expectations\new procedure_in progress\temp datasets\missing"

# for dataset in dataset_names:
#     for year in years:
#         df_name = f"{dataset}_{year}"
#         file_path = os.path.join(save_path, f"{df_name}.parquet")
#         if not os.path.exists(file_path):
#             print(f" File non trovato: {df_name}")
#             continue
#         df = pd.read_parquet(file_path)
#         # normalize Date
#         df['Date'] = pd.to_datetime(df['Date']).dt.normalize()
#         df_missing = df[df['Date'].isna()]
#         if not df_missing.empty:
#             df_missing = df_missing.assign(Dataset=df_name)
#             output_file = os.path.join(missing_path, f"missing_{df_name}.csv")
#             df_missing.to_csv(output_file, index=False)
#             dropped_instr = df_missing['Instrument'].dropna().unique().tolist()
#             print(f"\n Salvato: {output_file}")
#             print(f" Instruments con date mancanti: {dropped_instr}")
#         else:
#             print(f"✅ Nessun 'Date' mancante per {df_name}")

# Import missing all dataframe
# files = glob(os.path.join(missing_path, "missing_*.csv"))

# merged_missing = pd.DataFrame()
# for fp in files:
#     df = pd.read_csv(fp)
#     if "Dataset" not in df.columns:
#         dataset_name = os.path.basename(fp).replace("missing_", "").replace(".csv", "")
#         df.insert(0, "Dataset", dataset_name)
    
#     merged_missing = pd.concat([merged_missing, df], ignore_index=True)

# # Dataset as the first column 
# cols = ["Dataset"] + [c for c in merged_missing.columns if c != "Dataset"]
# merged_missing = merged_missing[cols]

# # save one comprehensive file reporting the missing values for all the datasets 
# output_file = os.path.join(missing_path, "missing_ALL_DATASETS.csv")
# merged_missing.to_csv(output_file, index=False)

##### Drop tickers that have more than 50 observations missing across the different datasets
##### Dropped tickers: 8303.T; 8355.T; 8604.T; ARLG.DE; CSGN.S; FRC.N; HBNK.CY; IL0A.I; PBCT.OQ^D22; SBNY.OQ; SIVB.OQ;  

# Save the cleaned datasets in cleaned_datasets folder - TO ELIMINATE
# out_path  = r"P:\ECB business areas\DGMF\SRF\Data management\Market data\IBES expectations\new procedure_in progress\cleaned_datasets"
# os.makedirs(out_path, exist_ok=True)

# Dataset base names
dataset_names = [
    "df_NetProfit", "df_ROE", "df_FeesandCommissionsInc",
    "df_Loans", "df_LoanLossProvisions", "df_NetIntInc",
    "df_TradingInc", "df_OperatingExp", "df_PreTaxProfit",
    "df_OPR", "df_TotalAssets", "df_TotalInc",
    "df_ShareholdersEquity", "df_TotalnonIntRev", "df_EPS",
    "df_DPS", "df_DPSYield"
    ]

tickers_to_drop = {"8303.T","8355.T","8604.T","ARLG.DE","CSGN.S","FRC.N",
                   "HBNK.CY","IL0A.I","PBCT.OQ^D22","SBNY.OQ","SIVB.OQ"}

for base in dataset_names:
    for year in years:
        name = f"{base}_{year}"
        df = globals()[name]  # get dataframe
        if "Instrument" in df.columns:
            globals()[name] = df[~df["Instrument"].isin(tickers_to_drop)]


# # check that unwanted tickers are gone from all dfs
# for base in dataset_names:
#     for year in years:
#         name = f"{base}_{year}"
#         if name not in globals():
#             print(f"{name} → skipped (not defined)")
#             continue
#         df = globals()[name]
#         if "Instrument" not in df.columns:
#             print(f"{name} → skipped (no Ticker column)")
#             continue
#         remaining = set(df["Instrument"]) & tickers_to_drop
#         if remaining:
#             print(f"{name} → ❌ FAILED, still found: {sorted(remaining)}")
#         else:
#             print(f"{name} → ✅ check passed, unwanted tickers removed")

### Change all the dataframes in long format for merging later
def wide_to_long (
    df,
    id_candidates=("Instrument", "Date", "Financial Period Absolute")):
    df = df.copy()
    # normalize Date if present
    if "Date" in df.columns:
        df["Date"] = pd.to_datetime(df["Date"], errors="coerce").dt.normalize()
    # split IDs vs values
    id_cols = [c for c in id_candidates if c in df.columns]
    value_cols = [c for c in df.columns if c not in id_cols]
    # melt
    long = df.melt(id_vars=id_cols, value_vars=value_cols,
                   var_name="Variable", value_name="Value")
    # split "Base - Statistic"
    extracted = long["Variable"].astype(str).str.extract(
        r"^(?P<Base>.*)\s*[-–—]\s*(?P<Statistic>[^-–—]+)$")
    long["Base"]      = extracted["Base"].fillna(long["Variable"]).str.strip()
    long["Statistic"] = extracted["Statistic"].str.strip()
    # now pivot so that Base becomes the column (e.g. Dividend Per Share)
    wide = long.pivot_table(
        index=id_cols + ["Statistic"],
        columns="Base",
        values="Value",
        aggfunc="first"
    ).reset_index()
    wide.columns.name = None
    return wide

for base in dataset_names:
    for yr in years:
        name = f"{base}_{yr}"
        if name in globals():
            try:
                globals()[f"{name}_stat"] = wide_to_long(globals()[name])
                print(f"{name} → ok")
            except Exception as e:
                print(f"{name} → error: {e}")
        else:
            print(f"{name} → skipped")

######## REMEMBER TO CHECK Total Non-Interest Revenue and Pre-Tax with merging 
# (2) each dataset is transformed into a long-format dataset to facilitate data visualization
# cleaned_path = r"P:\ECB business areas\DGMF\SRF\Data management\Market data\IBES expectations\new procedure_in progress\cleaned_datasets"
# long_path    = os.path.join(cleaned_path, "long_format")
# os.makedirs(long_path, exist_ok=True)
# files = glob(os.path.join(cleaned_path, "*.parquet"))
# for sample_fp in files:
#     print(f"\n Processing {os.path.basename(sample_fp)} ...")
#     # upload all the files
#     df = pd.read_parquet(sample_fp)
#     if "Date" in df.columns:
#         df["Date"] = pd.to_datetime(df["Date"], errors="coerce").dt.normalize()
#     # fix Pre-Tax → Pre Tax
#     df.columns = df.columns.str.replace("Pre-Tax", "Pre Tax", regex=False)
#     # fix Total Non-Interest Revenue → Total Non Interest Revenue
#     df.columns = df.columns.str.replace("Total Non-Interest Revenue", "Total Non Interest Revenue", regex=False)
#     # we identify ID
#     id_candidates = ["Instrument", "Date", "Financial Period Absolute"]
#     id_cols = [c for c in id_candidates if c in df.columns]
#     # all the others are variables
#     value_cols = [c for c in df.columns if c not in id_cols]
#     if not value_cols:
#         continue
#     # wide → long
#     df_long = df.melt(
#         id_vars=id_cols,
#         value_vars=value_cols,
#         var_name="Variable",
#         value_name="Value"
#     )
#     # separate Variable and Statistic
#     df_long[["Variable", "Statistic"]] = df_long["Variable"].str.rsplit(" - ", n=1, expand=True)
#     df_long = df_long[id_cols + ["Variable", "Statistic", "Value"]]
#     # save long format
#     out_long = os.path.join(long_path, os.path.basename(sample_fp).replace(".parquet", "_long.parquet"))
#     df_long.to_parquet(out_long, index=False)
#     # create wide format
#     df_wide = df_long.pivot_table(
#         index=id_cols + ["Statistic"],
#         columns="Variable",
#         values="Value"
#     ).reset_index()
#     df_wide.columns.name = None
#     # save wide format
#     out_wide = os.path.join(long_path, os.path.basename(sample_fp).replace(".parquet", "_wide.parquet"))
#     df_wide.to_parquet(out_wide, index=False)
#     print(f"✅ Salvato wide format: {out_wide}")

### Fix ROE statistic value
# Standard deviation variable report values in a duplicate NA column
def fix_roe_columns(df: pd.DataFrame) -> pd.DataFrame:
    """Merge 'Return on Equity' into 'Return On Equity' and drop the duplicate."""
    df = df.copy()
    low, proper = "Return on Equity", "Return On Equity"
    if low in df.columns and proper in df.columns:
        # Fill missing values in the proper-cased column from the other one
        df[proper] = df[proper].where(df[proper].notna(), df[low])
        df.drop(columns=[low], inplace=True)
    elif low in df.columns and proper not in df.columns:
        # Only the lower-cased version exists → rename to the proper-cased name
        df.rename(columns={low: proper}, inplace=True)
    return df
# Apply only to the ROE stat tables
for yr in years:
    name = f"df_ROE_{yr}_stat"
    if name in globals():
        try:
            globals()[name] = fix_roe_columns(globals()[name])
            print(f"{name} → ok")
        except Exception as e:
            print(f"{name} → error: {e}")
    else:
        print(f"{name} → skipped")

##### Merge all the dataframes for each projection year
def merge_all_for_year_left(year, dataset_names):
    # ROE is the anchor
    base_name = f"df_ROE_{year}_stat"
    if base_name not in globals():
        raise ValueError(f"{base_name} not found")
    base = globals()[base_name].copy()
    
    id_keys = ["Instrument", "Date", "Financial Period Absolute", "Statistic"]

    for name in dataset_names:
        if name == "df_ROE":
            continue  # already base
        stat_name = f"{name}_{year}_stat"
        if stat_name in globals() and isinstance(globals()[stat_name], pd.DataFrame):
            base = pd.merge(base, globals()[stat_name], on=id_keys, how="left")
        else:
            print(f"{stat_name} → skipped (not found)")
    
    # add a column for the year (FY1, FY2, FY3)
    base["Financial Period Relative"] = year
    return base

# final projection year dfs
FY1 = merge_all_for_year_left("FY1", dataset_names)
FY2 = merge_all_for_year_left("FY2", dataset_names)
FY3 = merge_all_for_year_left("FY3", dataset_names)

# # post-merge verification checks
# def check_merge(df, name):
#     print(f"\n--- Checking {name} ---")
    
#     # 1. Duplicates on merge keys
#     dups = df.duplicated(subset=["Instrument", "Date", "Financial Period Absolute", "Statistic"]).sum()
#     print(f"Duplicate rows on merge keys: {dups}")
    
#     # 2. Columns created by name collisions
#     xy_cols = [c for c in df.columns if c.endswith(("_x", "_y"))]
#     print(f"Columns with _x/_y: {xy_cols}")
    
#     # 3. Null values
#     nulls = df.isna().sum()
#     print("Columns with missing values (top 10):")
#     print(nulls[nulls > 0].sort_values(ascending=False).head(10))
    
#     print(f"Rows: {len(df)}, Cols: {df.shape[1]}")
#     print("--- Done ---")

# check_merge(FY1, "FY1")
# check_merge(FY2, "FY2")
# check_merge(FY3, "FY3")

### Create the final dataframe by concatenation
all_df = pd.concat([FY1, FY2, FY3], ignore_index=True)

# rename variables in line with suba harmonized 
rename_map = {
    "Return On Equity": "roe",
    "Net Income": "inc_net",
    "Fees & Commissions Income": "inc_fee_com_pas_oth",
    "Loans": "loans_net",
    "Loan Loss Provisions": "loans_exp_prov",
    "Net Interest Income": "inc_int_net",
    "Trading Income": "inc_trading",
    "Operating Expense": "exp_operating",
    "Operating Profit":"pro_operating",
    "Pre-Tax Profit":"pro_beforetax",
    "Total Assets": "tassets",
    "Total Income": "tincome",
    "Shareholders Equity": "equity",
    "Total Non-Interest Revenue": "no_int_income",
    "Dividend Per Share": "DPS",
    "Earnings Per Share": "EPS",
}

all_df = all_df.rename(columns=rename_map)

### Compute some new derived variables

# Other_expenses = tincome - (exp_operating + pro_operating)
mask = all_df[["tincome", "exp_operating", "pro_operating"]].notna().all(axis=1)
all_df["Other_expenses"] = (all_df["tincome"] - (all_df["exp_operating"] + all_df["pro_operating"])).where(mask, np.nan)

# Other_net_income = pro_beforetax + loans_exp_prov - pro_operating
mask = all_df[["pro_beforetax", "loans_exp_prov", "pro_operating"]].notna().all(axis=1)
all_df["Other_net_income"] = (all_df["pro_beforetax"] + all_df["loans_exp_prov"] - all_df["pro_operating"]).where(mask, np.nan)

# Taxes = pro_beforetax - inc_net
mask = all_df[["pro_beforetax", "inc_net"]].notna().all(axis=1)
all_df["Taxes"] = (all_df["pro_beforetax"] - all_df["inc_net"]).where(mask, np.nan)

# Other_operating_income = tincome - (inc_int_net + inc_fee_com_pas_oth + inc_trading)
mask = all_df[["tincome", "inc_int_net", "inc_fee_com_pas_oth", "inc_trading"]].notna().all(axis=1)
all_df["Other_operating_income"] = (
    all_df["tincome"] - (all_df["inc_int_net"] + all_df["inc_fee_com_pas_oth"] + all_df["inc_trading"])
).where(mask, np.nan)

###############################################################################
# Export final dataset
# devo.to_parquet(all_df, path="/dlb_dgmf/share/PROJECTS/DVP/outlets/fsr/banks/profitability_projections_refinitiv_ibes_ea.parquet")
all_df.to_parquet("P:/ECB business areas/DGMF/SRF/Data management/Market data/IBES expectations/profitability_projections_refinitiv_ibes_banks.parquet")
print(all_df.columns.tolist())