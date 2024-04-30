from dataclasses import dataclass
from geopy.distance import distance as geopy_distance

@dataclass
class Location():
    lat: float 
    lon: float 
        
    def __str__(self):
        return f"({self.lat:.2f},{self.lon:.2f})"
    
    def distance(self, other):
        """Calculate the distance between this location and another location."""
        coords1 = (self.lat, self.lon)
        coords2 = (other.lat, other.lon)
        dist_km = geopy_distance(coords1, coords2).km
        return dist_km
        
    def pretty_print(self):
        if self.lat < 0:
            lat_hemisphere = "S"
        else:
            lat_hemisphere = "N"
        
        if self.lon < 0:
            lon_hemisphere = "W"
        else:
            lon_hemisphere = "O"
        
        return f"{abs(self.lat):.2f}°{lat_hemisphere}, {abs(self.lon):.2f}°{lon_hemisphere}"  
        

