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
    
    def __str__(self):
        return "({:e} , {:e} , {:e})".format(self.x, self.y, self.z)

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
    compute_velocity(bodies, time_step )
    update_location(bodies, time_step)

def plot_output(bodies, outfile = None):
    fig = plot.figure()
    ax = fig.add_subplot(1,1,1, projection='3d')
    max_range = 0
    i=0 #indice pour les couleurs
    for current_body in bodies: 
        randcolor = (random.random(),random.random(), random.random())
        max_dim = max(max(current_body["x"]),max(current_body["y"]),max(current_body["z"]))
        if max_dim > max_range:
            max_range = max_dim
        ax.plot(current_body["x"], current_body["y"], current_body["z"], ':', c = randcolor, label = current_body["name"])        
        ax.plot(current_body["x"][-1],current_body["y"][-1], current_body["z"][-1], 'o', c = randcolor) 
        i+=1
    
    ax.set_xlim([-max_range,max_range])    
    ax.set_ylim([-max_range,max_range])
    ax.set_zlim([-max_range,max_range])
    ax.legend(title='legende', bbox_to_anchor=(1, 1), loc='upper left')


    if outfile:
        plot.savefig(outfile, dpi=300)
    else:
        plot.show()

# Calcule la distance minimal entre 2 astres
def dist(a,b):
    return math.sqrt(math.pow(a.x-b.x,2)+math.pow(a.y-b.y,2)+math.pow(a.z-b.z,2))

# Calcule une matrice de distances entre tout les astres d'une liste
def matriceDistance(bodies):
    matrice = []
    for body1 in bodies:
        line = []
        for body2 in bodies:
            line.append(dist(body1.location,body2.location))
        matrice.append(line)
    return matrice

# Affiche en notation scientifique toutes les distances entre des astres
def prettyDistance(bodies):
    matrice = matriceDistance(bodies)
    for index1,body1 in enumerate(bodies):
        print( "Distances retaives a l'astre " + body1.name + ':' )
        for index2,body2 in enumerate(bodies):
            if index1 != index2:
                print("     "+body2.name+" : {:e}".format(matrice[index1][index2]))
        print("\n")
            
# retourne un tuple x,y des indexes de la distance la plus faible dans une matrice
def minDistMatrice(matrice):
    currentMin = matrice[0][1]
    tupleIndexes = (0,1)
    for x in range(len(matrice)):
        line = matrice[x]
        for y in range(len(line)):
            if x!=y and matrice[x][y] < currentMin:
                currentMin = matrice[x][y]
                tupleIndexes = (x,y)
    return tupleIndexes

def run_simulation(bodies, names = None, time_step = 1, number_of_steps = 10000, report_freq = 100):
    #liste pour chaque corps
    body_locations_hist = []
    for current_body in bodies:
        body_locations_hist.append({"x":[], "y":[], "z":[], "name":current_body.name})
    
    minMatrice = matriceDistance(bodies)
    minSimTuple = minDistMatrice(minMatrice)
    minSim = minMatrice[minSimTuple[0]][minSimTuple[1]]

    for i in range(1,number_of_steps):
        compute_gravity_step(bodies, time_step)            
        
        matrice = matriceDistance(bodies)
        currMinTuple = minDistMatrice(matriceDistance(bodies))
        currMin = matrice[currMinTuple[0]][currMinTuple[1]]
        if currMin < minSim:
            minSim = currMin
            minSimTuple = currMinTuple
            minMatrice = matrice

        if i % report_freq == 0:
    

            for index, body_location in enumerate(body_locations_hist):
                body_location["x"].append(bodies[index].location.x)
                body_location["y"].append(bodies[index].location.y)           
                body_location["z"].append(bodies[index].location.z)       
    
    print("La distance minimal de la simulation est est {:e}".format(minMatrice[minSimTuple[0]][minSimTuple[1]])+" il sagit de la distance entre "+bodies[minSimTuple[0]].name+" et "+bodies[minSimTuple[1]].name )

    return body_locations_hist        
            
#planet data (location (m), mass (kg), velocity (m/s)
bodiesSpec = [
    { "id": 0, "name" : "Soleil" , "mass":2e30},
    { "id": 1, "name" : "Mercimam" , "mass":3.285e23},
    { "id": 2, "name" : "Ventdanus" , "mass":4.8e24},
    { "id": 3, "name" : "Kerbin" , "mass":6e24},
    { "id": 4, "name" : "MarsAtak" , "mass":2.4e24},
    { "id": 5, "name" : "Jupiter-3" , "mass":1e28},
    { "id": 6, "name" : "Satourne" , "mass":5.7e26},
    { "id": 7, "name" : "UrAnus" , "mass":8.7e25},
    { "id": 8, "name" : "Padidé :(" , "mass":1e26},
    { "id": 9, "name" : "PlutoLeIench" , "mass":1.3e22},
    { "id" : "apophis", "name" : "Armagedon" , "mass":5e10}
]

locationConstant = 1.5e11
velocityConstant = 1731460

if __name__ == "__main__":
    
    
    sim_start_date = "2021-01-01"
    time = Time(sim_start_date).jd
    
    bodies = []

    for bodySpec in bodiesSpec:
       
        obj = Horizons(id=bodySpec["id"], location="@sun", epochs=time, id_type='id').vectors()
        bodies.append(
            body( 
                location=point(
                    np.double(obj['x']*locationConstant),
                    np.double(obj['y']*locationConstant),
                    np.double(obj['z']*locationConstant)),
                mass = bodySpec["mass"],
                velocity = point(
                    np.double(obj['vx']*velocityConstant),
                    np.double(obj['vy']*velocityConstant),
                    np.double(obj['vz']*velocityConstant)),
                name = bodySpec["name"]
            )
        )

   
    print(prettyDistance(bodies))

    matrice = matriceDistance(bodies)
    minDist = minDistMatrice(matrice)

    print(bodies[0].velocity)

    print("La distance minimal est {:e}".format(matrice[minDist[0]][minDist[1]])+" il sagit de la distance entre "+bodies[minDist[0]].name+" et "+bodies[minDist[1]].name )

    motions = run_simulation(bodies, time_step = 1000, number_of_steps = 10000, report_freq = 100)
    
    plot_output(motions, outfile = 'orbits.png')