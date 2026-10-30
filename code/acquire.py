"""Retrieve public environmental snapshots; never modifies original data."""
from pathlib import Path
import urllib.request,urllib.parse,json,hashlib,datetime,zipfile,io
from concurrent.futures import ThreadPoolExecutor
import pandas as pd

ROOT=Path(__file__).resolve().parents[1]
AUDIT=ROOT.parent/'Shark_Storm_Overlap_Audit_2026-09-16'
SOURCE=ROOT/'sources';SOURCE.mkdir(exist_ok=True)
RESULT=ROOT/'results';RESULT.mkdir(exist_ok=True)

def fetch(url,name,limit=40_000_000):
    p=SOURCE/name
    if p.exists():return p
    row={'url':url,'file':name,'retrieved_utc':datetime.datetime.now(datetime.timezone.utc).isoformat()}
    try:
        req=urllib.request.Request(url,headers={'User-Agent':'Scientific reproducibility audit; public data'})
        with urllib.request.urlopen(req,timeout=40) as response:
            b=response.read(limit+1)
            if len(b)>limit:raise ValueError('Size limit exceeded')
            row.update(status=response.status,bytes=len(b),sha256=hashlib.sha256(b).hexdigest(),content_type=response.headers.get('Content-Type'),final_url=response.url)
        p.write_bytes(b)
    except Exception as e:row.update(status='failed',error=str(e))
    # Separate per-request log files avoid parallel append races.
    (SOURCE/(name+'.retrieval.json')).write_text(json.dumps(row,indent=2))
    print(name,row['status'],row.get('bytes',row.get('error')),flush=True)
    return p if p.exists() else None

def main():
    p=pd.read_csv(AUDIT/'results/screened_positions.csv.gz')
    p=p[p.deployment.eq('14_07-132980')].copy();p['time']=pd.to_datetime(p.time,utc=True)
    p=p[p.time.between('2015-05-01','2015-06-21')]
    p.to_csv(RESULT/'shark_modelled_positions.csv',index=False)
    print('Position bounds',p[['lat','lon']].agg(['min','max']).to_dict(),flush=True)
    # Explicit nominal grid locations, rounded to nearest ERA5 quarter degree.
    cells=p[['lat','lon']].mul(4).round().div(4).drop_duplicates().sort_values(['lat','lon'])
    cells.columns=['request_lat','request_lon'];cells['role']='trajectory'
    # Fixed-location scenarios around the event position, not statistical error bounds.
    scenarios=pd.DataFrame([{'request_lat':26.875+a,'request_lon':-113.15+b,'role':'spatial_scenario'} for a,b in [(0,0),(.5,0),(-.5,0),(0,.5),(0,-.5)]])
    cells=pd.concat([cells,scenarios],ignore_index=True);cells['site_id']=['site_%03d'%i for i in range(len(cells))]
    cells.to_csv(RESULT/'weather_requested_sites.csv',index=False)
    q=urllib.parse.urlencode({'latitude':','.join(cells.request_lat.astype(str)),'longitude':','.join(cells.request_lon.astype(str)),'start_date':'2015-05-01','end_date':'2015-06-20','hourly':'wind_speed_10m,pressure_msl,temperature_2m','wind_speed_unit':'ms','models':'era5','timezone':'GMT','cell_selection':'sea'})
    lo=(int((p.lat.min()-.75)*2)//1)/2;hi=(int((p.lat.max()+.75)*2)+1)/2
    left=(int((p.lon.min()-.75)*2)-1)/2;right=(int((p.lon.max()+.75)*2)+1)/2
    sl=f'[(2015-05-01T00:00:00Z):1:(2015-06-20T23:00:00Z)][({lo}):1:({hi})][({left}):1:({right})]'
    wave='https://erddap.aoml.noaa.gov/hdb/erddap/griddap/WaveWatch_2015.csv?'+urllib.parse.quote(','.join(v+sl for v in ['Thgt','Tper','shgt','whgt']),safe=',:()')
    jobs=[('https://archive-api.open-meteo.com/v1/archive?'+q,'era5_sites.json'),(wave,'ww3_region.csv'),('https://erddap.aoml.noaa.gov/hdb/erddap/griddap/WaveWatch_2015.das','ww3_metadata.txt'),('https://www.nhc.noaa.gov/data/tcr/EP022015_Blanca.pdf','NHC_Blanca.pdf')]
    with ThreadPoolExecutor(max_workers=4) as pool:list(pool.map(lambda a:fetch(*a),jobs))
    zpath=AUDIT/'sources/tags/14_07-132980.zip'
    with zipfile.ZipFile(zpath) as z:
        names=['out-RTC.csv','out-Summary.csv','sources.txt','out-Labels.csv']
        for n in names:
            b=z.read(n);(SOURCE/n).write_bytes(b)
        nc=next(n for n in z.namelist() if n.endswith('.nc'))
        (SOURCE/nc).write_bytes(z.read(nc))
    print('Acquisition finished',flush=True)

if __name__=='__main__':main()
