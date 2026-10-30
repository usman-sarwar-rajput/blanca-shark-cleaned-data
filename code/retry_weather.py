from acquire import SOURCE,RESULT,fetch
import pandas as pd,urllib.parse,json
from concurrent.futures import ThreadPoolExecutor

def main():
    sites=pd.read_csv(RESULT/'weather_requested_sites.csv')
    def get(r):
        q=urllib.parse.urlencode({'latitude':r['request_lat'],'longitude':r['request_lon'],'start_date':'2015-05-01','end_date':'2015-06-20','hourly':'wind_speed_10m,pressure_msl,temperature_2m','wind_speed_unit':'ms','models':'era5','timezone':'GMT','cell_selection':'sea'})
        return fetch('https://archive-api.open-meteo.com/v1/archive?'+q,r['site_id']+'.json')
    with ThreadPoolExecutor(max_workers=3) as pool:paths=list(pool.map(get,sites.to_dict('records')))
    print('Successful sites',sum(p is not None for p in paths),'of',len(paths))

if __name__=='__main__':main()
