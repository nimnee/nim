"""
http://sciience.tumblr.com/post/101026151217/quick-python-script-for-making-voronoi-polygons
"""
import os, osgeo, math, numpy as np
from osgeo import gdal, ogr, osr
from scipy.spatial import Delaunay
import scipy
root = 'E:/google_places/' # directory with input point shp
input_shp_full = root+'cafe.shp' # point shp
fieldUID = 'Field1' # field in point shp with pt IDs
out_dir = root+'out-voronoi/' # output directory to hold buffered points and Voronoi polygons
proj = 4326 #http://spatialreference.org/ref/epsg/wgs-84/

#################
### INPUT DATA ###
#################
print "fieldUID:", fieldUID
# read in shp file
print "- reading",input_shp_full
drv = ogr.GetDriverByName('ESRI Shapefile')
ptShp = drv.Open(input_shp_full)
ptLayer = ptShp.GetLayer(0)
print ptLayer
ptSRS = ptLayer.GetSpatialRef()
x_min, x_max, y_min, y_max = ptLayer.GetExtent()

# collect point coordinates and ID
ptList = []
ptDict = {}
for pt in ptLayer:
    #print "pt----------------",pt
    ID_index = pt.GetFieldIndex(fieldUID)
    ptID = pt.GetField(ID_index)
    #print "ptID----------------",ptID
    ptGeom = pt.GetGeometryRef()
    ptX = float(str(ptGeom).split(' ')[1].strip('('))
    ptY = float(str(ptGeom).split(' ')[2].strip(')'))
    ptDict[ptID] = [ptX,ptY] # a bit redundant to have dict and list, but list is used for Delaunay input
    ptList.append([ptX,ptY]) 

print "\t>> read",len(ptList),"points"
numPtList = np.array(ptList) # set-up for input to Delaunay   

#########################
### RADIAL BUFFER PTS ###
########################

# declare buffering distance in units appropriate for spatial reference system
# here, we're in WGS84 so we'll buffer by 1/100th of a degree
bufferDistance = 0.01 # degrees
print "- buffering",input_shp_full,"at",bufferDistance,"deg"

# set-up output buffered point shp
drv = ogr.GetDriverByName('Esri Shapefile')
buffShp = out_dir+input_shp_full.split('/')[-1].replace('.shp','_buffer.shp')
if os.path.exists(buffShp): os.remove(buffShp)
ds = drv.CreateDataSource(buffShp)

#srs = osr.SpatialReference()
#srs.ImportFromEPSG(4326)


layer = ds.CreateLayer('', None, ogr.wkbPolygon)#None,
layer.CreateField(ogr.FieldDefn('Id', ogr.OFTInteger))
defn = layer.GetLayerDefn()

# run through points and buffer each by established distance
ptCounter = 0
for each in ptDict:
    ptLon = ptDict[each][0]
    ptLat = ptDict[each][1]
    #print "buffering point at lat-lon",ptLat,"-",ptLon
    pt_wkt = "POINT ("+str(ptLon)+' '+str(ptLat)+')'
    pt = ogr.CreateGeometryFromWkt(pt_wkt)
    
    geom = pt.Buffer(bufferDistance)
    feat = ogr.Feature(defn)
    feat.SetField('Id', each)
    feat.SetGeometry(geom)  
    layer.CreateFeature(feat)
    feat = geom = None
    ptCounter+=1

print "\t>> buffered",ptCounter,'points in',
layer = ds = None

########################
### VORONOI POLYGONS ##
########################

# References
    # http://en.wikipedia.org/wiki/Circumscribed_circle#Circumscribed_circles_of_triangles
    # https://stackoverflow.com/questions/12374781/how-to-find-all-neighbors-of-a-given-point-in-a-delaunay-triangulation-using-sci
    # https://stackoverflow.com/questions/10650645/python-calculate-voronoi-tesselation-from-scipys-delaunay-triangulation-in-3d

# read in points to Delaunay function
print "- converting points in",input_shp_full,"to Voronoi polygons"
tri = Delaunay(numPtList)

# returns point locations of each triangle vertex
# tri.vertices: each row represents one simplex (triangle) in the triangulation,
# with values referencing indices of the input point list
p = tri.points[tri.vertices]
#print "p:",p
print "p:",len(p)
print "p:",p[0]
# find facets containing each input point
print '\tfinding Delaunay triangle facets'
triDict={}
i = 0
while i > 'exported '+str(ptCounter)+' Voronoi vertices':
    nodeLayer = outNodeShp = None
    print i

for i in xrange(len(p)):
    triDict[i]=p[i]
#print "triDict:",triDict
#print "ptDict:",ptDict
# build polygons with vertices at identified nodes
# assign input shp pt IDs to polygons
#vorIdDict = triDict
import collections


vorIdDict = collections.defaultdict(list)
for pt in triDict:
    #print "pt: ",pt
    meanLon0 = -1
    meanLat0 = -1
    if len(triDict[pt]) != 3:
        print "!=3"
    else:
        #meanLon0 = sum(node[0] for node in triDict[pt]) / len(triDict[pt])
        print triDict[pt][0]

        # meanLon0 = -1
        # meanLat0 = -1
        A = np.array(triDict[pt][0])
        B = np.array(triDict[pt][1])
        C = np.array(triDict[pt][2])
        a = np.linalg.norm(C - B)
        b = np.linalg.norm(C - A)
        c = np.linalg.norm(B - A)
        s = (a + b) / 2
        R = a * b / 4 / np.sqrt(s * (s - a) * (s - b) * (s - c))
        b1 = a * a * (b * b + c * c - a * a)
        b2 = b * b * (a * a + c * c - b * b)
        b3 = c * c * (a * a + b * b - c * c)
        P = np.column_stack((A, B, C)).dot(np.hstack((b1, b2, b3)))
        P /= b1 + b2 + b3

    for node in triDict[pt]:
        exist = False
        returnID = -1
        for ptID, ptArray in ptDict.iteritems():
            if ptArray[0] == node[0] and ptArray[1] == node[1]:
                exist = True
                returnID = ptID
        if exist:
            #vorIdDict[returnID].append([meanLon0, meanLat0])
            vorIdDict[returnID].append(P)
        else:
            print "Error! 142"

print "vorIdDict:",type(vorIdDict)
# https://stackoverflow.com/questions/1709283/how-can-i-sort-a-coordinate-list-for-a-rectangle-counterclockwise
# https://gamedev.stackexchange.com/questions/13229/sorting-array-of-points-in-clockwise-order
# https://en.wikipedia.org/wiki/Graham_scan
def sortCCW(node):
    return math.atan2(node[1] - meanLat, node[0] - meanLon)

# set-up output polygon shp
vorShp = out_dir+input_shp_full.split('/')[-1].replace('.shp','_voronoi.shp')
print "- creating output polygon shp",vorShp.split('/')[-1]
if os.path.exists(vorShp): os.remove(vorShp)
drv = ogr.GetDriverByName('ESRI Shapefile')
outShp = drv.CreateDataSource(vorShp)
layer = outShp.CreateLayer('', None,ogr.wkbPolygon)
layer.CreateField(ogr.FieldDefn('Id', ogr.OFTInteger))
layerDefn = layer.GetLayerDefn()

# find nodes surrounding polygon centroid
# sort nodes in counterclockwise order
# create polygon perimeter through nodes
print "- building Voronoi polygons around point..."
for pt in vorIdDict:
    #print "\t",pt
    #print "0"
    #print "type:",type(vorIdDict[pt][0])
    meanLon = sum(node[0] for node in vorIdDict[pt])/len(vorIdDict[pt])

    #print "1"
    meanLat = sum(node[1] for node in vorIdDict[pt])/len(vorIdDict[pt])
    #print "2"
    #hullList = vorIdDict[pt].tolist()
    hullList = vorIdDict[pt]
    #print "3"
    hullList.sort(key=sortCCW)
    #print "4"
  
    poly = ogr.Geometry(ogr.wkbPolygon)
    ring = ogr.Geometry(ogr.wkbLinearRing)
    i = 0
    for node in hullList:
        if i==0:
            loopLon = node[0] # grab first node to close ring
            loopLat = node[1]
        ring.AddPoint(node[0],node[1])
        i+=1
    ring.AddPoint(loopLon,loopLat)
    #ring.AddPoint(meanLon, meanLat)
    poly.AddGeometry(ring)
    feat = ogr.Feature(layerDefn)
    feat.SetField('Id', pt)
    feat.SetGeometry(poly)  
    layer.CreateFeature(feat)
    feat = poly = ring = None
layer = outShp = None