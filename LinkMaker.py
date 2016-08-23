# -*- coding: utf-8 -*-
#############################
# This code is designed to extract data from the 
# referenced Google spreadsheet and to convert that
# to a Shapefile of 'densified' lines that will 
# enable us to show great circle lines between the
# pair of institutions. 
#
# In the script, we are looking for is a way to take 
# a text label for location and to obtain coordinates 
# for it using Wikipedia.
#
# Once we have a coordinate pair (technically, *if* we 
# have a coordinate pair) for an institution and
# a partner then we try to draw a great circle
# route between the two. The line is broken into  
# multiline segments that ensure it won't get re-drawn
# as a straight line in a projection.
# 
# User-specified setting:
# -----------------------
# 1. Type of view
view = 'Default'  # Standard Robinson projection centered on the Atlantic
#view = 'Alternate' # Robinson-type projection centered on the Pacific (-lon0 = 150)
#############################

import os
import re
import math
import gspread 
import numpy as np
import requests
import wikipedia
from itertools import compress
import pandas as pd
import geopandas as gpd
from shapely.geometry import Point, LineString, MultiLineString
from geographiclib.geodesic import Geodesic

# Function to extract coordinates from
# Wikipedia web page. Uses wikipedia lib.
# To resolve a DisambiguationError you'll 
# have to look at the Wikipedia site in order
# to find out exactly what the best title is
# to ensure that one and only one page is 
# returned from the title parameter.
def getCoordinates(title):
	"""
	: title: Specify the unique Wikipedia title for the page to be retrieved
	"""
    # Figure it's helpful to say what we're looking for
    # in case things go wrong
    print("Fetching Wikipedia page for: " + title)
    
    # Create a Wikipedia parser
    pg = wikipedia.page(title)
    
    # Look for coordinates values in the parser
    try: 
        return pg.coordinates
    # Didn't find coordinate values
    except KeyError, e: 
        print("No coordinates found on page for: " + title)
        return None
    # There is more than one Wikipedia page matching this title
    except DisambiguationError, e: 
        print("Multiple pages with this title: " + title)
        return None

# Function to find Location URL on a
# Wikipedia page. Sometimes we don't have
# coordinates value, but we do have data 
# in the raw Wikipedia page that provides 
# location information. In this case we have
# to try to rip it out of the location URL
# embedded in the page.
def getLocation(title, pattern='locality'):
	"""
	: title: Specify the unique Wikipedia title for the page to be retrieved
	: pattern: defaults to "locality", but can be overridden if the coordinates are stored elsewhere
	"""
    # Figure it's helpful to say what we're doing
    # in case things go wrong
    print("Falling back to raw Wikipedia page: " + pg.url)
    
    # Create a Wikipedia parser
    pg = wikipedia.page(title)

    # Request the page and go looking for the
    # pattern using a regex-style search. Note
    # that we default to 'locality', which seems
    # to cover most eventualities.
    f      = requests.get(pg.url)
    coord  = f.text.find(pattern)
    # Extract the next 200 characters from the
    # starting point of our pattern match
    substr = f.text[(coord):(coord + 200)]
    
    # Now find the unique title for the page
    # that contains the locality information
    loc    = re.findall('title="([^"]+?)"', substr)
    
    # And request *that* so that we can grab
    # its coordinate data. We shouldn't have a 
    # Disambiguation error here since Wikipedia
    # shouldn't allow a locality page to be 
    # non-unique...
    pg = wikipedia.page(loc[0])
    try:
        return pg.coordinates
    except KeyError, e:
        print("Unable to find coordinates for " + loc[0])
        return None

# Function to search for best set of
# coordinates available -- they will
# either be associated with the object
# itself, or with the object referenced
# as a 'Location' somewhere on the page.
def findCoordinates(title):
	"""
	This should be called by most users as it is simply
	a wrapper function for the two _other_ ways of searching
	for coordinates on a Wikipedia page.
	: title: Specify the unique Wikipedia title for the page to be retrieved
	"""
    cPair = getCoordinates(title)
    if cPair is None:
        cLoc = getLocation(title)
        if cLoc is not None:
            cPair = cLoc
    return cPair

# How to calculate great circle lines on.
# Adapted from: http://geographiclib.sourceforge.net/html/python/examples.html#computing-waypoints
# Requires: GeographicLib

# Note that we pass in the coordinate X and Y as separate 
# parameters. This avoid confusion about whether we should
# have X then Y, or Y then X, depending on whether you're a
# geographer or a computer scientist.
def getGreatCircleFromFloats(starty, startx, endy, endx, geod=Geodesic.WGS84):
    """
    Function for calculating great circle routes. Uses
    floats for inputs; these are expected to be derived 
    from WGS84 coordinates because we've initialised geod
    as WGS84 Geodesic object.
    : starty: the starting latitude point in WGS84
    : startx: the starting longitude point in WGS84
    : endy: the ending latitude point in WGS84
    : endx: the ending longitude point in WGS84
    : geod: initialised by default as a Geodesic.WGS84 object
    """
    # If we couldn't get useable coordinates
    # then just return None since we can't 
    # actually draw the circle.
    if np.isnan(starty) or np.isnan(startx) or np.isnan(endy) or np.isnan(endx):
        return None
    
    # We don't particularly care about the Y, but having 
    # the startx be less than than that start y helps to 
    # ensure that the line gets drawn consistently. Otherwise
    # you might see most of the lines following a normal 
    # path and then one that runs much further to the East 
    # or West because it should be wrapping around and isn't.
    if (startx > endx):
        (tmpx, tmpy)     = (endx, endy)
        (endx, endy)     = (startx, starty)
        (startx, starty) = (tmpx, tmpy)
    
    # The bulk of the work is simply drawing this inverse
    # line now that we've got the coordinates sorted out.
    l = geod.InverseLine(float(starty), float(startx), float(endy), float(endx),
                    Geodesic.LATITUDE | Geodesic.LONGITUDE) # | Geodesic.LONG_UNROLL)
	
    # Initialise the arc information
    da = 1; n = int(math.ceil(l.a13 / da)); da = l.a13 / n

    # Try to deal with break at -180 (int'l dateline)
    ml = list() # list of line lists
    line = list() # empty list for a new line
    lastLon = None # track where we jump across the dateline

    for i in range(n + 1):
        a = da * i
        g = l.ArcPosition(a, Geodesic.LATITUDE |
                             Geodesic.LONGITUDE | Geodesic.LONG_UNROLL)

        # if view == 'Default':
        #     if lastLon is not None and ((g['lon2'] <= -180 and lastLon >= -180) or (g['lon2'] >= 180 and lastLon <= 180)):
        #         if len(line) > 1:
        #             #print("Breaking line at: " + str(line))
        #             ml.append(line)
        #             line = list()
        #             #print("--- break ---")
        # elif view == 'Alternate':
        #     if lastLon is not None and ((g['lon2'] <= -330 and lastLon >= -330) or (g['lon2'] >= 30 and lastLon <= 30)):
        #         if len(line) > 1:
        #             #print("Breaking line at: " + str(line))
        #             ml.append(line)
        #             line = list()
        #             #print("--- break ---")

        # Copy to a new variable to make it
        # easier to see what's happening
        x = g['lon2']
        y = g['lat2']
        # if view == 'Default': # Standard Robinson projection centered on Atlantic
        #     if x < -180:
        #         x = 180 - (abs(x) - 180)
        #     elif x > 180:
        #         x = -180 + (abs(x) - 180)
        # elif view == 'Alternate': # Robinson projection centered on Pacific
        #     if x < -330:
        #         x = 30 - (abs(x) - 330)
        #     elif x > 30:
        #         x = -330 + (abs(x) - 30)
        line.append(Point(x, y))
        lastLon = g['lon2']
        #print "{:.5f} {:.5f}".format(g['lat2'], g['lon2'])

    # Append the last line that was running
    # unless it's a single point, in which
    # case we want to stick it on to the
    # previous multi-line string
    if len(line) > 1:
        ml.append(line)
	
	# Now assemble this into a MultiLineString
	# object that can be written to a shapefile
    for i in range(len(ml)):
        ml[i] = LineString(ml[i])
        #LineString(lineElems)
    return MultiLineString(ml)

def getGreatCircleFromPoints(start, end, geod=Geodesic.WGS84):
    """
    Utility function that allows us to use Point 
    class objects instead of raw floats to calculate
    the great circle route.
    : start: a shapely.Point object
    : end: a shapely.Point object
    : geod: initialised by default as a Geodesic.WGS84 object
    """
    if type(start) is not Point:
        #print("Start is not a Point: " + str(type(start)))
        return None
    elif type(end) is not Point:
        #print("End is not a Point: " + str(type(end)))
        return None
    else:
        return getGreatCircleFromFloats(start.y, start.x, end.y, end.x, geod)

# I'm keeping this stuff in Dropbox for now
# -- may need to change
os.chdir(os.path.expanduser("~"))
if os.path.exists("./Dropbox/Global Health"):
    os.chdir("./Dropbox/Global Health")
else:
    os.chdir("./Documents/Dropbox/Global Health")

###############
# Load the appropriate XLSX from Dropbox
df = pd.read_excel('./Global Health Partnership Data - 20160607.xlsx')

# Reformat the column names to something
# that is easier to work with in Python
cols  = list(df.columns.values)
ncols = list()
for c in cols:
    c = c.strip()
    c = c.strip(',;?"!')
    ncols.append(re.sub(' ','_',re.sub('^partner ','p',c)))
#print(ncols)
df.columns = ncols
del(cols, ncols, c)

# Fill in the missing data by column
# using the forward fill feature of pandas.
# Notice that we can only do to this for
# columns that are meant to be fully
# populated.
for col in ['university','city','country','has_global_health_programme','has_global_health_partnership']:
    df[col] = df[col].ffill()

######
# Would be good to add a quick check to see if the
# Excel file has been modified since we last created
# the pickles below -- if so then we need to check
# that we've not got a meaningful change and update
# if necessary
######

# Now find the locations of all of the universities
# -- many will have this on the Wikipedia page, but
# a few of them list only the town/city.
pik1 = './Global Health Partnership Data - University Locations.pickle'
if os.path.isfile(pik1):
    print("Restoring university locations from pickle...")
    cdf = pd.read_pickle(pik1)
else:
    print("Finding university locations using Wikipedia...")
    locs = dict()
    unis = list(df.university.unique())
    for u in unis:
        if not locs.has_key(u):
            locs[u] = findCoordinates(u)

    # Create a data structure suitable for use as
    # a data frame.
    lats = list()
    lons = list()
    for key in unis:
        if len(locs[key]) > 1:
            lats.append(locs[key][0])
            lons.append(locs[key][1])
        else:
            lats.append(None)
            lons.append(None)

    d = {'university': unis, 'lat': lats, 'lon': lons}
    cdf = pd.DataFrame(data=d)
    cdf.to_pickle(pik1)

cdf.head()

pik2 = './Global Health Partnership Data - Partner Locations.pickle'
if os.path.isfile(pik2):
    print("Restoring partner locations from pickle...")
    pdf = pd.read_pickle(pik2)
else:
    print("Finding partner locations using Wikipedia...")
    locs = dict()

    partners = list(df.pcity.unique())
    for p in partners:
        if not locs.has_key(p) and isinstance(p, basestring):
            locs[p] = findCoordinates(p)

    # Create a data structure suitable for use as
    # a data frame.
    lats = list()
    lons = list()
    for p in partners:
        print(p)
        try:
            print(locs[p])
            if len(locs[p]) > 1:
                lats.append(locs[p][0])
                lons.append(locs[p][1])
            else:
                lats.append(None)
                lons.append(None)
        except KeyError, e: # Some NaNs kicking around
            pass

    # Check for NaNs
    toRemove   = list(pd.isnull(pd.Series(partners)))
    institutes = list(compress(partners, [not i for i in toRemove]))

    d = {'location': institutes, 'lat': lats, 'lon': lons}
    pdf = pd.DataFrame(data=d)
    pdf.to_pickle(pik2)

# And write out the partner df as a
# shapefile.
crs = None
geometry = [Point(xy) for xy in zip(pdf.lon, pdf.lat)]
pgdf     = gpd.GeoDataFrame(pdf, crs=crs, geometry=geometry)
pgdf.drop('lat', axis=1, inplace=True)
pgdf.drop('lon', axis=1, inplace=True)
pgdf.crs = {'init':'epsg:4326'} # Set up QGIS map in EPSG:3395?
pgdf.to_file('Global Health Partner Locations.shp')
pgdf.head()

# And write out the institutional df as a
# shapefile
geometry = [Point(xy) for xy in zip(cdf.lon, cdf.lat)]
cgdf     = gpd.GeoDataFrame(cdf, crs=crs, geometry=geometry)
cgdf.drop('lat', axis=1, inplace=True)
cgdf.drop('lon', axis=1, inplace=True)
cgdf.crs = {'init':'epsg:4326'} # Set up QGIS map in EPSG:3395?
cgdf.to_file('Global Health Institution Locations.shp')
cgdf.head()

gdft = pd.merge(df, cgdf, how='left', on='university')
gdft.head()

gdf = pd.merge(gdft, pgdf, how='left', left_on='pcity', right_on='location')
gdf.head()
del(gdft)

#getGreatCircleFromPoints(gdf.geometry_x[0], gdf.geometry_y[0]) # Should be an error if first institution has no partners
#getGreatCircleFromPoints(gdf.geometry_x[1], gdf.geometry_y[1]) # Should return a LineString as second institution has one or more partners

# Now get a linestring for every row in the data frame
def applyCircleRoute(row):
    #print('Creating great circle for ' + str(row.university) + ', ' + str(row.pinstitution))
    return getGreatCircleFromPoints(row.geometry_x, row.geometry_y)

if 'Path' in gdf.columns:
    gdf.drop('Path', axis=1, inplace=True)
gdf['Path'] = gdf.apply(applyCircleRoute, axis=1)
gdf.head()

#gdf.ix[gdf.university=='Yale University',]

# Now convert gdf to something that can be written
# out as a shapefile...
t = gdf.set_geometry('Path', crs={'init':'epsg:4326'})
t.drop('geometry_x', axis=1, inplace=True)
t.drop('geometry_y', axis=1, inplace=True)
t.head()
# Remove non-existent lines
t = t.ix[pd.notnull(t['pcity']),]
t.to_file('Global Health Great Circles ' + view + '.shp')
#with open('Global Health Great Circles.geojson', 'w') as f:
#    f.write(t.to_json())
