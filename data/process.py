import pandas as pd
from tqdm import tqdm
import dask.dataframe as dd
from dask.diagnostics import ProgressBar

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

occ_agg = [
    {'occ3_security' : ['occ3_guard', 'occ3_protect']},
    {'occ3_cleaners' : ['occ3_clean', 'occ3_janitor']},
    {'occ3_service' : ['occ3_food', 'occ3_shealth', 'occ3_recreation', 'occ3_child', 'occ3_others']},
    {'occ2_management' : ['occ2_exec', 'occ2_mgmtrel', 'occ2_cleric']},
    {'occ2_professional' : ['occ2_prof', 'occ2_tech']},
    {'occ2_sales' : ['occ2_finsales', 'occ2_retsales']},
    {'occ2_firepolice' : ['occ2_firepol']},
    {'occ2_primary' : ['occ2_farmer', 'occ2_otheragr', 'occ2_mining']},
    {'occ2_trades' : ['occ2_mechanic', 'occ2_constr', 'occ2_operator', 'occ2_transp', 'occ2_product', 'occ2_others']}
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

print('Processing data...')
occupation_mapping = {}
for category_dict in occ_agg:
    for category, occupations_list in category_dict.items():
        for occupation in occupations_list:
            occupation_mapping[occupation] = category

print('Reading data...')
df = dd.read_csv('raw/usa_00015.csv', assume_missing=True, dtype={'INDNAICS': 'string'})
met_codes = pd.read_csv('raw/met_codes.csv')
met_codes.set_index('code', inplace=True)
print('Mapping occpations...')
# df['METAREA'] = df['METAREA'].map(met_codes['METAREA'])
df['occupation'] = df['OCC1990'].apply(map_occupations)
# df['occupation'] = df['occupation'].map(occupation_mapping)

print('Mapping areas...')
df_cz = pd.read_csv('raw/cz_mappings.csv')
df_cz = df_cz[['LMA/CZ', 'FIPS']]
df['FIPS'] = df['STATEFIP'] * 1000 + df['COUNTYFIP']
df['COMZONE'] = df['FIPS'].map(df_cz.set_index('FIPS')['LMA/CZ'])

print('Mapping sectors...')
df['INDNAICS'] = df['INDNAICS'].astype(str).str[:2]
naics_codes = pd.read_csv('raw/2022_NAICS_codes.csv')
naics_codes = naics_codes[['Sector', 'Name']]
naics_codes.dropna(inplace=True)
naics_codes.set_index('Sector', inplace=True)
naics_codes.loc['0'] = 'N/A'
naics_codes.loc['99'] = 'Unemployed'
naics_codes.loc['50'] = 'Transportation and Warehousing'
naics_codes.loc['3M'] = 'Manufacturing'
naics_codes = expand_ranges(naics_codes)
df['INDNAICS'] = df['INDNAICS'].map(naics_codes['Name'])

years = df['YEAR'].unique().compute().astype(int)

total_pop = []
city_rent_all = pd.DataFrame()

area = 'COMZONE'

for year in tqdm(years, desc='Generating datasets: '):
    print(f'Processing data for year {year}...')
    print('Filtering data...')
    with ProgressBar():
        current = df[df['YEAR'] == year].compute(scheduler='threads')

    print('Processing...')
    current['INCWAGE'] = current['INCWAGE'] * current['HHWT']
    current['AVERAGE INCWAGE'] = current['INCWAGE'] / current['HHWT']
    current['RENT'] = current['RENT'] * current['HHWT']
    
    current['HH_RENT'] = current['RENT'].notnull().astype(int)
    current['HH_RENT'] = current['HH_RENT'] * current['HHWT']

    city_occ = current.pivot_table(index=area, columns='occupation', values='HHWT', aggfunc='sum')
    city_occ_wage = current.pivot_table(index=area, columns='occupation', values='AVERAGE INCWAGE', aggfunc='mean')
    city_sector_wage = current.pivot_table(index=area, columns='INDNAICS', values='INCWAGE', aggfunc='sum')
    sector_occ_wage = current.pivot_table(index=area, columns='occupation', values='INCWAGE', aggfunc='sum')

    city_rent = current.pivot_table(index=area, values='RENT', aggfunc='sum')
    if city_rent.empty:
        continue
    
    city_pop = current.pivot_table(index=area, values='HH_RENT', aggfunc='sum')
    city_rent[f'{year}'] = city_rent['RENT'] / city_pop['HH_RENT']

    city_rent_all = pd.concat([city_rent_all, city_rent[f'{year}']], axis=1)

    city_occ = city_occ[~city_occ.index.isin(['Not in identifiable area'])]
    try:
        city_sector_wage = city_sector_wage.drop(columns=['N/A', 'Unemployed'])
        sector_occ_wage = sector_occ_wage[~sector_occ_wage.index.isin(['N/A', 'Unemployed'])]
    except:
        pass

    print('Saving datasets...')
    city_occ.to_csv(f'processed/city_occ_employment/city_occ_e_{year}.csv')
    city_occ_wage.to_csv(f'processed/city_occ_wage/city_occ_w_{year}.csv')
    city_sector_wage.to_csv(f'processed/city_sec_wage/city_sec_w_{year}.csv')
    sector_occ_wage.to_csv(f'processed/sec_occ_wage/sec_occ_w_{year}.csv')
    # city_rent.to_csv(f'processed/city_rent/city_rent_{i}.csv')

    total_pop.append({'year': year, 'pop': current['HHWT'].sum()})

    print('Dumping cache...')
    del current, city_occ, city_occ_wage, city_sector_wage, city_rent, city_pop
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