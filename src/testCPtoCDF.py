import cPickle as cp
import BubbleModelStorage as BMS

print("Load CPickle File")
denList = cp.load(open('v2', 'r'))

print("Create a denList by Tempo")
denList = BMS.find_by_tempo(denList, -2.25010)

print("Creating NETCDF")
BMS.create_netcdf(denList, name="dataset_by_tempo")
