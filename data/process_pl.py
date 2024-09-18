import pandas as pd
from tqdm import tqdm
import polars as pl

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
df = pl.scan_csv('raw/usa_00015.csv')
met_codes = pd.read_csv('raw/met_codes.csv')
met_codes.set_index('code', inplace=True)

print('Mapping occpations...')
df = df.with_columns(
    pl.col('OCC1990').map_elements(map_occupations, return_dtype=pl.Utf8).alias('occupation')
)

print('Creating area mappings...')
df_cz = pd.read_csv('raw/cz_mappings.csv')
df_cz = df_cz[['LMA/CZ', 'FIPS']]
cz_names = pd.read_csv('../data/raw/archive/cz_county.csv')
cz_names = cz_names.dropna()
cz_names['LMA/CZ'] = cz_names['LMA/CZ'].astype(str).str[:-2].astype(int)
idx = cz_names.groupby(['LMA/CZ'])['Labor Force'].idxmax()
cz_names = cz_names.loc[idx].reset_index(drop=True)
cz_names = cz_names[['LMA/CZ', 'County Name']].set_index('LMA/CZ')
cz_names['County Name'] = cz_names['County Name'].str.replace('"', '')
df_cz = df_cz.merge(cz_names, left_on='LMA/CZ', right_index=True)
df = df.with_columns(
    (pl.col('STATEFIP').cast(pl.Int32) * 1000 + pl.col('COUNTYFIP').cast(pl.Int32)).alias('FIPS')
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
    current['COMZONE'] = current['FIPS'].map(df_cz.set_index('FIPS')['County Name'])

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
print('Done!')