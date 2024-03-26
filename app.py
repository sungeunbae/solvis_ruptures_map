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
import plotly.express as px
import plotly.graph_objects as go


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

all_sol = InversionSolution.from_archive(cru_archive)
fault_names = all_sol.fault_sections.ParentName.unique().tolist()


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
    rupt_ids = rupt_ids_within_rate_range(all_sol, rate_range[0], rate_range[1])
    print(rupt_ids)
    sol2 = InversionSolution().filter_solution(all_sol, rupt_ids)
    if location is not None:
        rupt_ids = within_radius_rupt_ids(sol2, location, radius) # radius is in meters. 100km->100e3
        sol3 = InversionSolution().filter_solution(sol2, rupt_ids)
    else:
        sol3 = sol2

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

def generate_folium_map(sol, rupt_ids, location, radius, fmap=None,rupt_id=0):


    all_ruptures_fg = folium.FeatureGroup(name="All ruptures", overlay=False,control=False,show=True)
    all_ruptures_df = section_participation(sol,rupt_ids) 
    draw_rupture(all_ruptures_df, all_ruptures_fg)
    all_ruptures_fg.add_to(fmap)

    ruptures_fg = []
    rupture_df=section_participation(sol, [rupt_ids[rupt_id]])
    fg=folium.FeatureGroup(name=f"Scenario {rupt_id+1}", overlay=True)
    draw_rupture(rupture_df, fg, "red")
    fg.add_to(fmap)

    if location is not None:
        folium.Circle(location=(cities[location][0],cities[location][1]),radius=radius).add_to(fmap) # add circle
    folium.LatLngPopup().add_to(fmap)
    folium.LayerControl(collapsed=False,draggable=True).add_to(fmap)

    return fmap

@st.cache_data
def convert_df(df):
    if df is not None:
        return df.to_csv().encode('utf-8')
    else:
        return ""


def plot_scatter(all_ruptures_df, filtered_df=None):
   
    # Create the scatter plot
    fig = go.Figure()

    # Add all ruptures data points (green markers)
    fig.add_trace(go.Scatter(x=all_ruptures_df["Magnitude"], y=all_ruptures_df["Annual Rate"],
                             mode='markers', name=f"All Ruptures ({len(all_ruptures_df)})",
                             marker=dict(color='green', size=10, opacity=0.7)))

    if filtered_df is not None:
    # Add filtered data points (red markers)
        fig.add_trace(go.Scatter(x=filtered_df["Magnitude"], y=filtered_df["Annual Rate"],
                                 mode='markers', name=f"Selected Ruptures ({len(filtered_df)})",
                                 marker=dict(color='red', size=8, opacity=0.5)))

        # Add ellipse shape
        fig.add_shape(type="circle",
                      xref="x", yref="y",
                        x0=min(filtered_df["Magnitude"]), y0=min(filtered_df["Annual Rate"]),
                        x1=max(filtered_df["Magnitude"]), y1=max(filtered_df["Annual Rate"]),
                        opacity=0.9,
                      line=dict(color="red"))

    # Set y-axis to be logarithmic
    fig.update_yaxes(title_text="Annual Rate", type="log")
    fig.update_xaxes(title_text="Magnitude")
    
    fig.update_layout(width=1200, height=800)

    # Show the interactive plot
    #fig.show()
    return fig



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

    st.session_state.fmap = folium.Map(location=[-42.1, 172.8], zoom_start=6, tiles='cartodbpositron')  # Centered around New Zealand
    # Input widgets
    faults_to_include = st.sidebar.multiselect("Faults (optional)", fault_names, placeholder="Faults (optional)",default=None)
    location = st.sidebar.selectbox("Locations", sorted(cities.keys()),index=None)

    if location:
        radius_enabled = True

    radius = st.sidebar.selectbox("Radius (km)", (10,20,30,40,50,100,200), disabled=not radius_enabled)
    magnitude_range = st.sidebar.slider("Magnitude", 6,10, (6,10))
    rate_range = st.sidebar.slider("Rate (1eN/yr)", -20, 0, (-20,0))

    all_ruptures = all_sol.ruptures_with_rupture_rates
    if "sol" in st.session_state:
        
        selected_ruptures = st.session_state.sol.ruptures_with_rupture_rates
    else:
        selected_ruptures = None

    fig = plot_scatter(all_ruptures, selected_ruptures)
    st.plotly_chart(fig)
    st.download_button(label="all_ruptures.csv",data=convert_df(all_ruptures), file_name='all_ruptures.csv', mime='text/csv')
    st.download_button(label="selected_ruptures.csv",data=convert_df(selected_ruptures), file_name='selected_ruptures.csv', mime='text/csv',disabled=(selected_ruptures is None))


    scenario_val = st.slider("Scenarios",min_value=min_scenario,max_value=max_scenario,key="scenario", disabled=not scenario_enabled)
    if scenario_val == 0:
        scenario_val = 1


    print(faults_to_include)
    print(location)
    print(radius)
    print(magnitude_range)
    print(rate_range)

    def call_get_ruptures():
        sol, rupt_ids = get_ruptures(faults_to_include, location, radius*1000, magnitude_range, (10**rate_range[0],10**rate_range[1]))
        st.session_state.max_scenario = len(rupt_ids)
        st.session_state.sol=sol
        st.session_state.rupt_ids=rupt_ids

        
    #    fig = plot_scatter(all_sol.ruptures_with_rupture_rates, sol.ruptures_with_rupture_rates)

    #    st.plotly_chart(fig)
#        all_sol.ruptures_with_rupture_rates.to_csv("all_ruptures.csv")
#        sol.ruptures_with_rupture_rates.to_csv("AF_ruptures.csv")

        

#    def call_generate_fmap():
#        st.session_state.fmap = generate_folium_map(st.session_state.sol, st.session_state.rupt_ids, location, radius*1000, fmap=st.session_state.fmap, rupt_id=scenario_val-1)  

    st.sidebar.button("Get ruptures", on_click=call_get_ruptures)

    
    if  st.button("Generate Map", key="generate_map",disabled = not scenario_enabled):
        try:
            st.session_state.fmap = generate_folium_map(st.session_state.sol, st.session_state.rupt_ids, location, radius*1000, fmap=st.session_state.fmap, rupt_id=scenario_val-1)
        except Exception as e:
            raise e
        else:
            rr=st.session_state.sol.ruptures_with_rupture_rates
            rupt_id = st.session_state.rupt_ids[scenario_val-1]
            rupture = rr[rr["Rupture Index"]==rupt_id]

            st.sidebar.write(f"Rupure {scenario_val} of {len(st.session_state.rupt_ids)}")
            st.sidebar.write(f"Mean Rate: {rupture['Annual Rate'].values[0]:.2e} per year")
            st.sidebar.write(f"Magnitude: {rupture['Magnitude'].values[0]}")
            st.sidebar.write(f"Area: {rupture['Area (m^2)'].values[0]/1e6} km2")
            st.sidebar.write(f"Length: {rupture['Length (m)'].values[0]/1e3} m")

    folium_static(st.session_state.fmap, width=1200, height=1200)
if __name__ == "__main__":
    main()
