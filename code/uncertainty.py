from acquire import ROOT,SOURCE,RESULT
import sys,json
sys.path.insert(0,str(ROOT/'dependencies'))
import h5py,numpy as np,pandas as pd

def main():
    with h5py.File(next(SOURCE.glob('*.nc'))) as f:
        lat=f['latitude'][:];lon=f['longitude'][:];times=pd.to_datetime(f['twelve_hour_timestamps'][:],unit='s',utc=True)
        rows=[];cells=[]
        for i,t in enumerate(times):
            if not pd.Timestamp('2015-05-01',tz='UTC')<=t<pd.Timestamp('2015-06-21',tz='UTC'):continue
            a=f['twelve_hour_likelihoods'][i].astype(float);a=np.where(a>=0,a,0);total=a.sum();assert total>0
            a/=total;order=np.argsort(a.ravel())[::-1];n=np.searchsorted(np.cumsum(a.ravel()[order]),.95)+1
            chosen=order[:n];yy,xx=np.unravel_index(chosen,a.shape)
            rows.append(dict(time=t,source_likelihood_sum=total,normalized_95_mass=a.ravel()[chosen].sum(),cells_95=n,lat_min=lat[yy].min(),lat_max=lat[yy].max(),lon_min=lon[xx].min(),lon_max=lon[xx].max(),southern_grid_boundary_mass=a[-1,:].sum()))
            # All nonzero model support, with normalized weight; these are model-conditional.
            for y,x in zip(*np.nonzero(a>0)):
                cells.append(dict(time=t,lat=float(lat[y]),lon=float(lon[x]),weight=float(a[y,x]),in_95=(y*a.shape[1]+x) in chosen))
    d=pd.DataFrame(rows);d.to_csv(RESULT/'location_uncertainty.csv',index=False)
    pd.DataFrame(cells).to_csv(RESULT/'location_likelihood_cells.csv.gz',index=False)
    event=d[d.time.between('2015-06-05','2015-06-11')]
    print(event.to_string(index=False),flush=True)

if __name__=='__main__':main()
