import requests
import astropy.units as u
from astropy.time import Time
from astropy.coordinates import SkyCoord, EarthLocation, AltAz
from timezonefinder import TimezoneFinder
import pytz
from datetime import datetime


def convert_local_to_utc(location, local_time_str):
    """
    Converts a local time string to UTC based on an EarthLocation.

    Args:
        location (EarthLocation): The astropy location object.
        local_time_str (str): Format 'YYYY-MM-DD HH:MM:SS'
    """
    # 1. Extract Lat/Lon from EarthLocation
    lat = location.lat.deg
    lon = location.lon.deg
    print("*********************here************")
    # 2. Find the timezone name at these coordinates
    tf = TimezoneFinder()
    tz_name = tf.timezone_at(lng=lon, lat=lat)

    if not tz_name:
        return "Could not determine timezone for these coordinates."

    # 3. Create a timezone-aware datetime object
    local_tz = pytz.timezone(tz_name)
    naive_dt = datetime.strptime(local_time_str, "%Y-%m-%d %H:%M:%S")
    local_dt = local_tz.localize(naive_dt)

    # 4. Convert to UTC using Astropy
    # We pass the timezone-aware datetime directly to Time
    print("local TIME")
    print(local_dt)
    t = Time(local_dt)
    print(local_dt)
    print(t)
 #---------------------------------------------

    # 1. Get current UTC time in Astropy
    t_utc = Time.now()

    # 2. Get your local timezone info
    local_tz = datetime.now().astimezone().tzinfo

    # 3. Convert to a timezone-aware local datetime
    t_local = t_utc.to_datetime(timezone=local_tz)

    print(f"UTC Time:   {t_utc}")
    print(t.utc)
    print(f"Local Time: {t_local}")
 #   local_time = add_timezone_by_abbr(utc_time, "PST")
#  return t_utc.utc
    return t.utc

def get_location_by_ip():
    try:
        # Use ip-api.com (free, no API key required for low volume)
        response = requests.get("http://ip-api.com/json/")
        #response = requests.get('https://api.ipify.org?format=json')
        data = response.json()

        if data['status'] == 'success':
            lat = data['lat']
            lon = data['lon']
            city = data['city']
            region = data['regionName']

            # Create the Astropy EarthLocation object
            # Note: height is set to 0 as IP geo usually doesn't provide altitude
            location = EarthLocation(lat=lat * u.deg, lon=lon * u.deg, height=0 * u.m)

            print(f"Detected Location: {city}, {region}")
            print(f"Coordinates: {lat}, {lon}")
            return location
        else:
            print("Failed to get location data.")
            return None

    except Exception as e:
        print(f"An error occurred: {e}")
        return None







def get_m31_position():
    # Generate the object
    my_location = get_location_by_ip()
    print(my_location.lat)
    print(my_location.lon)
    obj_target = input("Enter the catalog or name of target: ")
    obs_time = input("Enter observation DateTime(YYYY-MM-DD HH:MM:SS): ")
    utc_time = convert_local_to_utc(my_location, obs_time)

    # 1. Get the current time in UTC
    current_time = Time.now()
    print(current_time)

    # 2. Set your location (Latitude, Longitude, Elevation)
    # Replace these with your actual coordinates
    # Example: Los Angeles, CA
   # my_location = EarthLocation(lat=34.0522*u.deg, lon=-118.2437*u.deg, height=71*u.m)

    # 3. Get M31 RA and Dec (ICRS frame)
    # This queries Sesame/SIMBAD to get the fixed celestial coordinates
    m31 = SkyCoord.from_name(obj_target)

    # 4. Transform to AltAz (Horizontal) frame for your location/time
    altaz_frame = AltAz(obstime=utc_time, location=my_location)
    m31_altaz = m31.transform_to(altaz_frame)

    # 5. Print the results
    print(f"Current UTC Time: {current_time.iso}")
    print("-" * 30)
    print(f"Target Coordinates at Observation Time: {obj_target}")
    print(f"Azimuth:  {m31_altaz.az.deg:.2f}°")
    print(f"Altitude: {m31_altaz.alt.deg:.2f}°")
    print(f"RA:     {m31.ra.hms.h:02.0f}h {m31.ra.hms.m:02.0f}m {m31.ra.hms.s:05.2f}s")
    print(f"Dec:    {m31.dec.dms.d:02.0f}d {abs(m31.dec.dms.m):02.0f}m {abs(m31.dec.dms.s):05.2f}s")
    print("-" * 30)
    print(f"Azimuth:  {m31_altaz.az.deg:.2f}°")
    print(f"Altitude: {m31_altaz.alt.deg:.2f}°")

    a1 = m31_altaz.az
    a2 = m31_altaz.alt

#===================================================================
    print("Slew to RA/Dec for preview now.")
    #altaz_coord = SkyCoord(alt=a2, az=a1, frame=altaz_frame)
    altaz_coord = SkyCoord(az=a1, alt=a2,
                         frame='altaz',
                         obstime=current_time,
                         location=my_location)
    crd = altaz_coord.transform_to('icrs')  # Or 'fk5', 'galactic', etc.
    print(f"RA:     {crd.ra.hms.h:02.0f}h {crd.ra.hms.m:02.0f}m {crd.ra.hms.s:05.2f}s")
    print(f"Dec:    {crd.dec.dms.d:02.0f}d {abs(crd.dec.dms.m):02.0f}m {abs(crd.dec.dms.s):05.2f}s")

if __name__ == "__main__":
    get_m31_position()