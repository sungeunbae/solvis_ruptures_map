# map_app.py
from pathlib import Path
import geopandas as gpd
from solvis import InversionSolution,section_participation,rupt_ids_above_rate,circle_polygon,export_geojson
from solvis.inversion_solution.typing import InversionSolutionProtocol

import streamlit as st
import folium
from streamlit_folium import st_folium, folium_static

from shapely.geometry import Polygon,LineString, LinearRing
from matplotlib.cm import ScalarMappable
import numpy as np
from matplotlib import colors
import matplotlib.pyplot as plt

from qcore import geo



R_EARTH = 6378.139


#cmap = plt.get_cmap('inferno')
#norm = plt.Normalize(1e-7,1e-4)

from matplotlib.colors import ListedColormap, BoundaryNorm

my_colors = ['#000000', '#4F0D6C', '#781C6D', '#D24644', '#ED6925', '#F6D746', '#FCFFA4']
boundaries = [7.5e-8, 2.5e-7, 7.5e-7, 2.5e-6, 7.5e-6, 2.5e-5, 7.5e-5]
cmap = ListedColormap(my_colors)
norm = BoundaryNorm(boundaries, cmap.N, clip=True)



WORK_DIR = Path("NSHM_data")
cru_archive = WORK_DIR / "CRU_fault_system_solution.zip"

sol = InversionSolution.from_archive(cru_archive)
fault_names = sol.fault_sections.ParentName.unique().tolist()


cities={'Wellington': (-41.276825, 174.777969),
 'Gisborne': (-38.662334, 178.017654),
 'Christchurch': (-43.52565, 172.639847),
 'Invercargill': (-46.413056, 168.3475),
 'Dunedin': (-45.8740984, 170.5035755),
 'Napier': (-39.4902099, 176.917839),
 'New Plymouth': (-39.0579941, 174.0806474),
 'Palmerston North': (-40.356317, 175.6112388),
 'Nelson': (-41.2710849, 173.2836756),
 'Blenheim': (-41.5118691, 173.9545856),
 'Whakatane': (-37.9519223, 176.9945977),
 'Greymouth': (-42.4499469, 171.2079875),
 'Queenstown': (-45.03, 168.66),
 'Auckland': (-36.848461, 174.763336),
 'Rotorua': (-38.1446, 176.2378),
 'Taupo': (-38.6843, 176.0704),
 'Whangarei': (-35.7275, 174.3166),
 'Levin': (-40.6218, 175.2866),
 'Tauranga': (-37.687, 176.1654),
 'Timaru': (-44.3904, 171.2373),
 'Oamaru': (-45.0966, 170.9714),
 'Pukekohe': (-37.2004, 174.901),
 'Hamilton': (-37.7826, 175.2528),
 'Lower Hutt': (-41.2127, 174.8997)}


def within_radius_rupt_ids(sol, cityname, radius):
    lat = cities[cityname][0]
    lon = cities[cityname][1]
    polygon = circle_polygon(radius_m=radius, lat=lat, lon=lon)
    #print(polygon)
    rupts = sol.get_ruptures_intersecting(polygon)

    return rupts

def rupt_ids_within_rate_range(sol: InversionSolutionProtocol, lower: float, upper: float, rate_column: str = "Annual Rate"):
    rr = sol.rupture_rates
    return rr[(rr[rate_column] >= lower) & (rr[rate_column] < upper) ]["Rupture Index"].unique()


def get_fault_planes(traces, dbottom, dtop, dip, dip_dir):
    planes = [] # each element is a list of 4 coordinates (lat,lon)
    lon1, lat1 = traces[0]
    lon2, lat2 = traces[1]
    strike = geo.ll_bearing(lon1, lat1, lon2, lat2, midpoint=True)


    if 180 > dip_dir - strike >= 0:
        # If the dipdir is not to the right of the strike, turn the fault around
        indexes = range(len(traces))
    else:
        indexes = range(len(traces) - 1, -1, -1)

    for i, i2 in zip(indexes[:-1], indexes[1:]):
        lon1, lat1 = traces[i]
        lon2, lat2 = traces[i2]

        height = dbottom - dtop
        width = abs(height / np.tan(np.deg2rad(dip)))

        bot_lat1, bot_lon1 = geo.ll_shift(lat1, lon1, width, dip_dir)
        bot_lat2, bot_lon2 = geo.ll_shift(lat2, lon2, width, dip_dir)

        acc = 4

        planes.append(
            [
                [round(lat1, acc), round(lon1, acc)],
                [round(lat2, acc), round(lon2, acc)],
                [round(bot_lat2, acc), round(bot_lon2, acc)],
                [round(bot_lat1, acc), round(bot_lon1, acc)],
            ]
        )

    return planes


def draw_rupture(faults_df, layer, color=None):

    if color is None:
        
        line_colors = [colors.to_hex(cmap(norm(faults_df.iloc[i]["Annual Rate"]))) for i in range(len(faults_df))]
        plane_color = "gray"
    else:
        plane_color = color

    for i in range(len(faults_df)):
        f=faults_df.iloc[i]
        
        fault_lines= f.geometry
        lon, lat = fault_lines.xy #fault_lines is LineString
        traces=list(zip(lon.tolist(),lat.tolist()))
        planes= get_fault_planes(traces,f.LowDepth,f.UpDepth,f.DipDeg,f.DipDir)
 
        new_planes = []
        for plane in planes:
            new_plane = []
            for lat,lon in plane:
                new_plane.append([lon,lat])
    
            new_planes.append(new_plane)
    
        if len(new_planes) == 0:
            continue
    
        ring_coords=[]
        for j, plane in enumerate(new_planes):
            ring_coords=ring_coords[:j]+plane+ring_coords[-j:]
    
        ring_geom = Polygon(ring_coords)
        ring = gpd.GeoDataFrame(index=[f.FaultID], crs='epsg:4326', geometry=[ring_geom])
        folium.GeoJson(ring,style_function=lambda x: 
                       {'fillColor': plane_color,'color': plane_color, "fillOpacity": 0.5,"weight":1 }, 
                       highlight_function=lambda feature: 
                       {"fillcolor": "yellow", "color":"green"}, 
                       name=f['ParentName'],
                       tooltip=f"{f['FaultName']} ({int(f.FaultID)}: {f['Annual Rate']:.2e}) ").add_to(layer) #,


        if color is None:
            folium.GeoJson(fault_lines,style_function=lambda x, i=i: {"color": line_colors[i], "weight":3}).add_to(layer)
        else:
            folium.GeoJson(fault_lines,style_function=lambda x: {"color": color, "weight":3}).add_to(layer)
    
    return layer

def get_ruptures(faults_to_include, location, radius, magnitude_range, rate_range):
    rupt_ids = rupt_ids_within_rate_range(sol, rate_range[0], rate_range[1])
    print(rupt_ids)
    sol2 = InversionSolution().filter_solution(sol, rupt_ids)
    rupt_ids = within_radius_rupt_ids(sol2, location, radius) # radius is in meters. 100km->100e3
    sol3 = InversionSolution().filter_solution(sol2, rupt_ids)

    if len(faults_to_include) > 0:
        rupt_ids_to_include = []
        for fault in faults_to_include:
            rupt_ids_with_fault = sol3.get_ruptures_for_parent_fault(fault)
            rupt_ids_to_include.extend(rupt_ids_with_fault)
        rupt_ids_to_include=np.array(sorted(set(rupt_ids_to_include)))

        sol4 = InversionSolution().filter_solution(sol3, rupt_ids_to_include) 
    else:
        sol4 = sol3 # include all
        rr = sol4.rupture_rates
        rupt_ids_to_include = rr["Rupture Index"].unique()
        print(sol3) 

    print(rupt_ids_to_include)
    print(sol4.ruptures_with_rupture_rates)

    return sol4, rupt_ids_to_include

def generate_folium_map(faults_to_include, location, radius, magnitude_range, rate_range, fmap=None, sol=None, rupt_ids=None, rupt_id=0):

    if sol is None:
        sol, rupt_ids = get_ruptures(faults_to_include, location, radius, magnitude_range, rate_range)

    if fmap is None:
        fmap = folium.Map(location=[-42.1, 172.8], zoom_start=6, tiles='cartodbpositron')  # Centered around New Zealand
        all_ruptures_fg = folium.FeatureGroup(name="All ruptures", overlay=False,control=False,show=True)
        all_ruptures_df = section_participation(sol,rupt_ids) 
        draw_rupture(all_ruptures_df, all_ruptures_fg)
        all_ruptures_fg.add_to(fmap)

    ruptures_fg = []
    rupture_df=section_participation(sol, [rupt_ids[rupt_id]])
    fg=folium.FeatureGroup(name=f"Scenario {rupt_id}", overlay=True)
    draw_rupture(rupture_df, fg, "red")
    fg.add_to(fmap)


    folium.Circle(location=(cities[location][0],cities[location][1]),radius=radius).add_to(fmap) # add circle
    folium.LatLngPopup().add_to(fmap)
    folium.LayerControl(collapsed=False,draggable=True).add_to(fmap)

    return fmap, sol, rupt_ids




def main():
    # Display the map
    st.set_page_config(layout="wide")
    st.title("Rupture Explorer")

    radius_enabled = False
    scenario_enabled = False
    print(st.session_state)

    if "max_scenario" not in st.session_state:
        min_scenario = 0
        max_scenario = 1
    else:
        min_scenario = 1
        max_scenario = st.session_state.max_scenario
        scenario_enabled = True


    # Input widgets
    faults_to_include = st.sidebar.multiselect("Faults (optional)", fault_names, placeholder="Faults (optional)",default=None)
    location = st.sidebar.selectbox("Locations", sorted(cities.keys()),index=None)

    if location:
        radius_enabled = True

    radius = st.sidebar.selectbox("Radius (km)", (10,20,30,40,50,100,200), disabled=not radius_enabled)
    magnitude_range = st.sidebar.slider("Magnitude", 6,10, (6,10))
    rate_range = st.sidebar.slider("Rate (1eN/yr)", -20, 0, (-20,0))

    scenario = st.sidebar.slider("Scenarios",min_value=min_scenario,max_value=max_scenario,key="num_scenarios", disabled=not scenario_enabled)
    if scenario == 0:
        scenario = 1

    # Initialize scenario 

    print(faults_to_include)
    print(location)
    print(radius)
    print(magnitude_range)
    print(rate_range)

    # Process user input
    if st.sidebar.button("Generate Map"):
        try:
            # Convert locations input to a list of tuples (latitude, longitude, city name)


            fmap,sol,rupt_ids = generate_folium_map(faults_to_include, location, radius*1000, magnitude_range, (10**rate_range[0],10**rate_range[1]),rupt_id=scenario-1)


            #fmap = folium.Map(location=[39.949610, -75.150282], zoom_start=16)
    #        st_data = st_folium(fmap, width=1500, height=1500)
            if len(rupt_ids)>0:
                st.session_state.max_scenario = len(rupt_ids)

            folium_static(fmap, width=1000, height=1200)

            #do something to refresh
        except Exception as e:
            st.error(f"Error: {e}")

    st.sidebar.write("Selected scenario:", scenario)
if __name__ == "__main__":
    main()
