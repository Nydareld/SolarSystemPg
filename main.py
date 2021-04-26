import math
import numpy as np
import random
import matplotlib.pyplot as plot
from mpl_toolkits.mplot3d import Axes3D


from astropy.time import Time
from astroquery.jplhorizons import Horizons
# rajout

plot.rcParams['figure.dpi'] = 300

class point:
    def __init__(self, x,y,z):
        self.x = x
        self.y = y
        self.z = z

class body:
    def __init__(self, location, mass, velocity, name = ""):
        self.location = location
        self.mass = mass
        self.velocity = velocity
        self.name = name

def calculate_single_body_acceleration(bodies, body_index):
    G_const = 6.67408e-11 #m3 kg-1 s-2
    acceleration = point(0,0,0)
    target_body = bodies[body_index]
    for index, external_body in enumerate(bodies):
        if index != body_index:
            r = (target_body.location.x - external_body.location.x)**2 + (target_body.location.y - external_body.location.y)**2 + (target_body.location.z - external_body.location.z)**2
            r = math.sqrt(r)
            tmp = G_const * external_body.mass / r**3
            acceleration.x += tmp * (external_body.location.x - target_body.location.x)
            acceleration.y += tmp * (external_body.location.y - target_body.location.y)
            acceleration.z += tmp * (external_body.location.z - target_body.location.z)

    return acceleration

def compute_velocity(bodies, time_step = 1):
    for body_index, target_body in enumerate(bodies):
        acceleration = calculate_single_body_acceleration(bodies, body_index)

        target_body.velocity.x += acceleration.x * time_step
        target_body.velocity.y += acceleration.y * time_step
        target_body.velocity.z += acceleration.z * time_step 


def update_location(bodies, time_step = 1):
    for target_body in bodies:
        target_body.location.x += target_body.velocity.x * time_step
        target_body.location.y += target_body.velocity.y * time_step
        target_body.location.z += target_body.velocity.z * time_step

def compute_gravity_step(bodies, time_step = 1):
    compute_velocity(bodies, time_step = time_step)
    update_location(bodies, time_step = time_step)

def plot_output(bodies, outfile = None):
    fig = plot.figure()
    colours = ['r','b','g','y','m','c','w','r','b','y']
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_range = 0
    i=0 #indice pour les couleurs
    for current_body in bodies: 
        max_dim = max(max(current_body["x"]),max(current_body["y"]),max(current_body["z"]))
        if max_dim > max_range:
            max_range = max_dim
        ax.plot(current_body["x"], current_body["y"], current_body["z"], ':', c = colours[i], label = current_body["name"])        
        ax.plot(current_body["x"][-1],current_body["y"][-1], current_body["z"][-1], 'o', c = colours[i]) 
        i+=1
    
    ax.set_xlim([-max_range,max_range])    
    ax.set_ylim([-max_range,max_range])
    ax.set_zlim([-max_range,max_range])
    ax.legend()        

    if outfile:
        plot.savefig(outfile, dpi=300)
    else:
        plot.show()

def run_simulation(bodies, names = None, time_step = 1, number_of_steps = 10000, report_freq = 100):

    #liste pour chaque corps
    body_locations_hist = []
    for current_body in bodies:
        body_locations_hist.append({"x":[], "y":[], "z":[], "name":current_body.name})
        
    for i in range(1,number_of_steps):
        compute_gravity_step(bodies, time_step = 1000)            
        
        if i % report_freq == 0:
            for index, body_location in enumerate(body_locations_hist):
                body_location["x"].append(bodies[index].location.x)
                body_location["y"].append(bodies[index].location.y)           
                body_location["z"].append(bodies[index].location.z)       

    return body_locations_hist        
            
#planet data (location (m), mass (kg), velocity (m/s)
sun = {"location":point(0,0,0), "mass":2e30, "velocity":point(0,0,0)}
mercury = {"location":point(0,5.7e10,0), "mass":3.285e23, "velocity":point(47000,0,0)}
venus = {"location":point(0,1.1e11,0), "mass":4.8e24, "velocity":point(35000,0,0)}
earth = {"location":point(0,1.5e11,0), "mass":6e24, "velocity":point(30000,0,0)}
mars = {"location":point(0,2.2e11,0), "mass":2.4e24, "velocity":point(24000,0,0)}
jupiter = {"location":point(0,7.7e11,0), "mass":1e28, "velocity":point(13000,0,0)}
saturn = {"location":point(0,1.4e12,0), "mass":5.7e26, "velocity":point(9000,0,0)}
uranus = {"location":point(0,2.8e12,0), "mass":8.7e25, "velocity":point(6835,0,0)}
neptune = {"location":point(0,4.5e12,0), "mass":1e26, "velocity":point(5477,0,0)}
pluto = {"location":point(0,3.7e12,0), "mass":1.3e22, "velocity":point(4748,0,0)}

if __name__ == "__main__":
    
    
    sim_start_date = "2021-01-01"
    time = Time(sim_start_date).jd
    

    obj1 = Horizons(id=1, location="@sun", epochs=time, id_type='id').vectors()
    obj2 = Horizons(id=2, location="@sun", epochs=time, id_type='id').vectors()
    obj3 = Horizons(id=3, location="@sun", epochs=time, id_type='id').vectors()
    obj4 = Horizons(id=4, location="@sun", epochs=time, id_type='id').vectors()
    obj5 = Horizons(id=5, location="@sun", epochs=time, id_type='id').vectors()
    obj6 = Horizons(id=6, location="@sun", epochs=time, id_type='id').vectors()
    obj7 = Horizons(id=7, location="@sun", epochs=time, id_type='id').vectors()
    obj8 = Horizons(id=8, location="@sun", epochs=time, id_type='id').vectors()
    obj9 = Horizons(id='apophis', location="@sun", epochs=time, id_type='id').vectors()

    #liste planète
    bodies = [
        body( location = sun["location"], mass = sun["mass"], velocity = sun["velocity"], name = "Soleil"),

        body( location=point(np.double(obj1['x']*1.5e11), np.double(obj1['y']*1.5e11),np.double(obj1['z']*1.5e11)), mass = mercury["mass"], velocity = point(np.double(obj1['vx']*1731460), np.double(obj1['vy']*1731460), np.double(obj1['vz']*1731460)), name = "Mercure"),
        body( location=point(np.double(obj2['x']*1.5e11), np.double(obj2['y']*1.5e11),np.double(obj2['z']*1.5e11)), mass = venus["mass"], velocity = point(np.double(obj2['vx']*1731460), np.double(obj2['vy']*1731460), np.double(obj2['vz']*1731460)), name = "Vénus"),
        body( location=point(np.double(obj3['x']*1.5e11), np.double(obj3['y']*1.5e11),np.double(obj3['z']*1.5e11)), mass = earth["mass"], velocity = point(np.double(obj3['vx']*1731460), np.double(obj3['vy']*1731460), np.double(obj3['vz']*1731460)), name = "Terre"),
        body( location=point(np.double(obj4['x']*1.5e11), np.double(obj4['y']*1.5e11),np.double(obj4['z']*1.5e11)), mass = mars["mass"], velocity = point(np.double(obj4['vx']*1731460), np.double(obj4['vy']*1731460), np.double(obj4['vz']*1731460)), name = "Mars"),
        #body( location=point(np.double(obj5['x']*1.5e11), np.double(obj5['y']*1.5e11),np.double(obj5['z']*1.5e11)), mass = jupiter["mass"], velocity = point(np.double(obj5['vx']*1731460), np.double(obj5['vy']*1731460), np.double(obj5['vz']*1731460)), name = "Jupiter"),
        #body( location=point(np.double(obj6['x']*1.5e11), np.double(obj6['y']*1.5e11),np.double(obj6['z']*1.5e11)), mass = saturn["mass"], velocity = point(np.double(obj6['vx']*1731460), np.double(obj6['vy']*1731460), np.double(obj6['vz']*1731460)), name = "Saturne"),
        #body( location=point(np.double(obj7['x']*1.5e11), np.double(obj7['y']*1.5e11),np.double(obj7['z']*1.5e11)), mass = uranus["mass"], velocity = point(np.double(obj7['vx']*1731460), np.double(obj7['vy']*1731460), np.double(obj7['vz']*1731460)), name = "uranus"),
        #body( location=point(np.double(obj8['x']*1.5e11), np.double(obj8['y']*1.5e11),np.double(obj8['z']*1.5e11)), mass = neptune["mass"], velocity = point(np.double(obj8['vx']*1731460), np.double(obj8['vy']*1731460), np.double(obj8['vz']*1731460)), name = "neptune"),
        body( location=point(np.double(obj9['x']*1.5e11), np.double(obj9['y']*1.5e11),np.double(obj9['z']*1.5e11)), mass = 5e10, velocity = point(np.double(obj9['vx']*1731460), np.double(obj9['vy']*1731460), np.double(obj9['vz']*1731460)), name = "Apophis"),
        ]

    
    
    motions = run_simulation(bodies, time_step = 100, number_of_steps = 8000, report_freq = 100)
    plot_output(motions, outfile = 'orbits.png')