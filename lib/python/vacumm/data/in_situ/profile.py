class Profile() :
    
    glider_begin = "EXGL"
    
    glider_type = "GLIDER"
    recopeca_type = "RECOPESCA"
    
    #cle des differentes tables
    time_key = "time"
    lon_key = "lon"
    lat_key = "lat"
    depth_key = "depth"
    temp_key = "temp"
    sal_key = "sal"


    """ Profile vertical """
    def __init__(self, time, platform_name, lon, lat):   
        
        self.name = "Vertical Profile"
        self.shortname = "VertProf"   
        
        self.time = time
        self.lon = lon
        self.lat = lat
        
        #tableau numpy
        self.depth = []
        self.temp = []
        self.sal = []
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...
                
        self.platform_name = platform_name
        if self.platform_name.startswith(self.glider_begin):
            self.platform_type = self.glider_type
        else:
            self.platform_type = self.recopeca_type


    # initialisation : creation des list du profil
    def add_time_step(self, depth=None, temp=None, sal=None):
        self.depth.append(depth)
        self.temp.append(temp)
        self.sal.append(sal)
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...

    
    #fin initialistion :creation table numpy
    def convert_tables(self):
        from numpy import array
        self.depth = array(self.depth, dtype='float')
        self.temp = array(self.temp, dtype='float')
        self.sal = array(self.sal, dtype='float')
        #...
        # AJOUTER NOUVELLES VARIABLES ICI
        #...

    
    # creation de tables pour tracer les profils
    # retourne une map de ndarray        
    def create_tabs(self):
        from numpy import ones_like
        from matplotlib.dates import date2num
        
        #duplication sur le nombre de profondeur (utile pour tracer le profil) 
        time_array = ones_like(self.depth) * date2num(self.time)
        lon_array = ones_like(self.depth) * float(self.lon)
        lat_array = ones_like(self.depth) * float(self.lat)        
        
        map_return = {self.time_key : time_array,
                      self.lon_key :lon_array,
                      self.lat_key :lat_array
                      }
        
        return map_return

    
    # retourne les infos date,lon et lat du profil
    def get_position(self):
        return [self.time, self.lon, self.lat]

