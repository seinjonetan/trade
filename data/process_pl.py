import pandas as pd
from tqdm import tqdm
import polars as pl
import glob
import os
import re

mappings = [
    ((405, 408), 'occ3_clean'),
    ((415, 427), 'occ3_protect'),
    ((415, 415), 'occ3_guard'),
    ((425, 427), 'occ3_guard'),
    ((433, 444), 'occ3_food'),
    ((445, 447), 'occ3_shealth'),
    ((448, 455), 'occ3_janitor'),
    ((457, 458), 'occ3_beauty'),
    ((459, 467), 'occ3_recreation'),
    ((468, 468), 'occ3_child'),
    ((469, 472), 'occ3_others'),
    ((3, 22), 'occ2_exec'),
    ((23, 37), 'occ2_mgmtrel'),
    ((43, 200), 'occ2_prof'),
    ((203, 235), 'occ2_tech'),
    ((243, 258), 'occ2_finsales'),
    ((274, 283), 'occ2_retsales'),
    ((303, 389), 'occ2_cleric'),
    ((417, 423), 'occ2_firepol'),
    ((473, 475), 'occ2_farmer'),
    ((479, 498), 'occ2_otheragr'),
    ((503, 549), 'occ2_mechanic'),
    ((558, 599), 'occ2_constr'),
    ((614, 617), 'occ2_mining'),
    ((628, 699), 'occ2_product'),
    ((703, 799), 'occ2_operator'),
    ((803, 889), 'occ2_transp'),
]

def map_occupations(occ1990dd):
    for (start, end), occupation in mappings:
        if start <= occ1990dd <= end:
            return occupation
    return 'other'

def expand_ranges(df):
    expanded_rows = []
    for index, row in df.iterrows():
        if '-' in index:
            start, end = map(int, index.split('-'))
            for i in range(start, end + 1):
                expanded_rows.append((str(i), row['Name']))
        else:
            expanded_rows.append((index, row['Name']))
    return pd.DataFrame(expanded_rows, columns=['Sector', 'Name']).set_index('Sector')

print('Reading data...')
df = pl.scan_csv('raw/usa_00019.csv')
# met_codes = pd.read_csv('raw/met_codes.csv')
# met_codes.set_index('code', inplace=True)

print('Mapping occpations...')
df = df.with_columns(
    pl.col('OCC1990').map_elements(map_occupations, return_dtype=pl.Utf8).alias('occupation')
)

print('Creating area mappings...')
# df_cz = pd.read_csv('raw/cz_mappings.csv')
# df_cz = df_cz[['LMA/CZ', 'FIPS']]
df_cz_1990 = pd.read_csv('raw/cw_puma1990_czone.csv', encoding='latin1')
df_cz_2000 = pd.read_csv('raw/cw_puma2000_czone.csv', encoding='latin1')
df_cz_1980 = pd.read_stata('raw/cw_ctygrp1980_czone_corr.dta')
# cz_names = pd.read_csv('../data/raw/archive/cz_county.csv')
# cz_names = cz_names.dropna()
# cz_names['LMA/CZ'] = cz_names['LMA/CZ'].astype(str).str[:-2].astype(int)
# idx = cz_names.groupby(['LMA/CZ'])['Labor Force'].idxmax()
# cz_names = cz_names.loc[idx].reset_index(drop=True)
# cz_names = cz_names[['LMA/CZ', 'County Name']].set_index('LMA/CZ')
# cz_names['County Name'] = cz_names['County Name'].str.replace('"', '')
# df_cz = df_cz.merge(cz_names, left_on='LMA/CZ', right_index=True)
df = df.with_columns(
    (pl.col('STATEFIP').cast(pl.Int32) * 10000 + pl.col('PUMA').cast(pl.Int32)).alias('FIPS'),
    (pl.col('STATEFIP').cast(pl.Int32) * 1000 + pl.col('CNTYGP98').cast(pl.Int32)).alias('ctygrp1980'),
)

print('Creating sector mappings...')
df = df.with_columns(
    pl.col('INDNAICS').cast(pl.Utf8).str.slice(0, 2).alias('INDNAICS')
)
naics_codes = pd.read_csv('raw/2022_NAICS_codes.csv')
naics_codes = naics_codes[['Sector', 'Name']]
naics_codes.dropna(inplace=True)
naics_codes.set_index('Sector', inplace=True)
naics_codes.loc['0'] = 'N/A'
naics_codes.loc['99'] = 'Unemployed'
naics_codes.loc['50'] = 'Transportation and Warehousing'
naics_codes.loc['3M'] = 'Manufacturing'
naics_codes = expand_ranges(naics_codes)

years = df.filter(pl.col('YEAR').is_not_null()).select(pl.col('YEAR').unique()).collect().to_series()
# years = pd.Series(years)[years >= 2010].to_list()

total_pop = []
city_rent_all = pd.DataFrame()

area = 'COMZONE'

for year in tqdm(years, desc='Generating datasets: '):
    print(f'Processing data for year {year}...')
    print('Filtering data...')
    current = df.filter(pl.col('YEAR') == year).collect().to_pandas()

    current['INCWAGE'] = current['INCWAGE'] * current['HHWT']
    current['AVERAGE INCWAGE'] = current['INCWAGE'] / current['HHWT']
    current['RENT'] = current['RENT'] * current['HHWT']

    current['HH_RENT'] = current['RENT'].notnull().astype(int)
    current['HH_RENT'] = current['HH_RENT'] * current['HHWT']

    print('Applying maps...')
    current['INDNAICS'] = current['INDNAICS'].map(naics_codes['Name'])
    # current['COMZONE'] = current['FIPS'].map(df_cz.set_index('FIPS')['County Name'])
    if year == 1990:
        current = current.merge(df_cz_1990, left_on='FIPS', right_on='puma1990', how='left')
        current['COMZONE'] = current['czone']
    elif year == 1980:
        current = current.merge(df_cz_1980, on='ctygrp1980', how='left')
        current['COMZONE'] = current['czone']
    else:
        current = current.merge(df_cz_2000, left_on='FIPS', right_on='puma2000', how='left')
        current['COMZONE'] = current['czone']

    city_occ = current.pivot_table(index=area, columns='occupation', values='HHWT', aggfunc='sum')
    city_occ_wage = current.pivot_table(index=area, columns='occupation', values='AVERAGE INCWAGE', aggfunc='mean')
    city_sector = current.pivot_table(index=area, columns='INDNAICS', values='HHWT', aggfunc='sum')
    sector_occ_wage = current.pivot_table(index='INDNAICS', columns='occupation', values='INCWAGE', aggfunc='sum')
    city_occ_wb = current.pivot_table(index=area, columns='occupation', values='AVERAGE INCWAGE', aggfunc='sum')

    city_rent = current.pivot_table(index=area, values='RENT', aggfunc='sum')
    if city_rent.empty:
        continue
    
    city_pop = current.pivot_table(index=area, values='HH_RENT', aggfunc='sum')
    city_rent[f'{year}'] = city_rent['RENT'] / city_pop['HH_RENT']

    if city_rent_all.empty:
        city_rent_all = city_rent[f'{year}']
        city_rent_all = pd.DataFrame(city_rent_all)
    else:
        city_rent_all = city_rent_all.merge(city_rent[f'{year}'], left_index=True, right_index=True, how='outer')

    city_occ = city_occ[~city_occ.index.isin(['Not in identifiable area'])]
    try:
        city_sector = city_sector.drop(columns=['N/A', 'Unemployed'])
        sector_occ_wage = sector_occ_wage[~sector_occ_wage.index.isin(['N/A', 'Unemployed'])]
    except:
        pass

    print('Saving datasets...')
    city_occ.to_csv(f'processed/city_occ_employment/city_occ_e_{year}.csv')
    city_occ_wage.to_csv(f'processed/city_occ_wage/city_occ_w_{year}.csv')
    city_sector.to_csv(f'processed/city_sec_employment/city_sec_e_{year}.csv')
    sector_occ_wage.to_csv(f'processed/sec_occ_wage/sec_occ_w_{year}.csv')
    city_occ_wb.to_csv(f'processed/city_occ_wb/city_occ_wb_{year}.csv')

    total_pop.append({'year': year, 'pop': current['HHWT'].sum()})

    print('Dumping cache...')
    del current, city_occ, city_occ_wage, city_sector, city_rent, city_pop

total_pop = pd.DataFrame(total_pop)
total_pop.to_csv('processed/total_pop.csv', index=False)
city_rent_all.to_csv('processed/city_rent.csv')
print('Success!')

print('Processing TFP data...')
tfp = pd.read_excel('raw/tfp.xlsx', header=1, skiprows=1)
tfp = tfp[tfp['Measure'] == 'Total factor productivity']
tfp = tfp[tfp['Units'].str.contains('Index')]
tfp['two_digit_naics'] = tfp['NAICS'].str[:2]
tfp.drop(columns=['Industry', 'NAICS', 'Basis', 'Measure', 'Units'], inplace=True)
tfp = tfp.groupby('two_digit_naics').mean()
tfp.rename(index={'MN': '33', 'DM': '33', 'ND': '33'}, inplace=True)
tfp.index = tfp.index.map(naics_codes['Name'])
tfp = tfp.groupby(tfp.index, axis=0).mean()
tfp.to_csv('processed/tfp.csv')

print('Creating master dataset...')
def get_data(directory, field_name):
    files = glob.glob(directory)
    files = [f for f in files if re.search(r'_(1980|199[0-9]|20[0-9]{2})\.csv$', f)]

    data = pd.DataFrame()

    for file in files:
        year = int(os.path.basename(file)[-8:-4])
        current = pd.read_csv(file)
        current['city_total'] = current.iloc[:, 1:].sum(axis=1)
        current = current.melt(id_vars=['COMZONE'], var_name='Occupation', value_name=field_name)
        current['Year'] = year
        data = pd.concat([data, current], ignore_index=True)

    return data

df_cw = pd.read_csv('raw/cw_puma2000_czone.csv', encoding='latin1')
df_names = pd.read_csv('raw/puma_names.csv')
df_names['puma2000'] = df_names['State10'] * 10000 + df_names['PUMA10']
df_names = df_names[['puma2000', 'PUMA10_Name']]
df_cw = df_cw.merge(df_names, on='puma2000', how='left')
df_cw = df_cw.sort_values(by=['czone', 'afactor'])
df_cw = df_cw.drop_duplicates(subset=['czone'])
df_cw = df_cw[['czone', 'PUMA10_Name']]
df_cw = df_cw.rename(columns={'czone': 'COMZONE', 'PUMA10_Name': 'NAME'})

employment = get_data('processed/city_occ_employment/*.csv', 'Employed')
employment = employment.merge(df_cw, on='COMZONE', how='left')
employment = employment[['Year', 'COMZONE', 'NAME', 'Occupation', 'Employed']]

wage = get_data('processed/city_occ_wage/*.csv', 'Wage')
wb = get_data('processed/city_occ_wb/*.csv', 'Wage_Bill')

final = employment.merge(wage, on=['Year', 'COMZONE', 'Occupation'], how='left')
final = final.merge(wb, on=['Year', 'COMZONE', 'Occupation'], how='left')
final.to_csv('master.csv', index=False)

print('Done!')