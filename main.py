
from vpython import *
from matplotlib.pyplot import *

G = 6.674e-11 # Newton gravitational constant
dt=1
earth_mass=5.9722e24
grav_param=(earth_mass+1e4)*G

class Orbit:
    def __init__(self,param,e,fi=0):
        self.p=param
        self.e=e
        self.fi=fi
        if self.e !=0:
            self.a=self.p/(1-self.e**2)
            self.b=sqrt(self.a*self.p)
        else:
            self.a=self.p
            self.b=self.p
    def draw_orbit(self):
        t=np.linspace(0,360,360)
        x=[(self.p/(1+self.e*cos(radians(i))))*cos(radians(i)) for i in t]
        y=[(self.p/(1+self.e*cos(radians(i))))*sin(radians(i)) for i in t]
        plot(x,y)
        show()


class Spacecraft:
    def __init__(self, mass, orbit, fi=0):
        self.mass = mass
        self.orbit = orbit
        self.r = np.array([orbit.p * cos(radians(fi)) / (1 + orbit.e * cos(radians(fi))),
                           orbit.p * sin(radians(fi)) / (1 + orbit.e * cos(radians(fi))), 0])
        self.v_n = sqrt(grav_param / orbit.p) * (1 + orbit.e * cos(radians(fi)))
        self.v_r = sqrt(grav_param / orbit.p) * orbit.e * sin(radians(fi))
        self.velocity = np.array([self.v_r * cos(radians(fi)) + self.v_n * cos(radians(fi + 90)),
                                  self.v_r * sin(radians(fi)) + self.v_n * sin(radians(fi + 90)), 0])
        self.h = np.linalg.norm(self.velocity) ** 2 - 2 * grav_param / np.linalg.norm(self.r)

    def simulation(self):
        scene.forward = vector(0, -1, -1)
        earth = sphere(pos=vector(0, 0, 0), radius=6.4e6, color=color.green)
        earth.mass = 6e24
        spacecraft = sphere(pos=vector(self.r[0], self.r[1], self.r[2]), radius=5e5, color=color.red,
                            make_trail=True, interval=10, retain=1000)

        dt = 10
        while True:
            rate(200)
            F = -G * earth.mass * self.mass * self.r / np.linalg.norm(self.r) ** 3
            self.velocity = self.velocity + F * dt / self.mass
            self.r = self.velocity * dt + self.r
            spacecraft.pos = vector(self.r[0], self.r[1], self.r[2])

    def get_momentum(self, delta):
        self.velocity = delta / self.mass + self.velocity

    def goman(self, new_orbit):
        pass


class Interval:
    def __init__(self, start_time, duration, spacecraft):
        self.interval = []
        self.speeds = []
        self.start_r = spacecraft.r
        self.start_time = start_time
        self.end_time = start_time + duration
        self.start_momentum = spacecraft.velocity * spacecraft.mass
        self.start_velocity = spacecraft.velocity
        spacecraft_p = spacecraft.velocity * spacecraft.mass
        spacecraft_pos = spacecraft.r
        for i in range(int(duration / dt)):
            F = -G * earth_mass * spacecraft.mass * spacecraft_pos / np.linalg.norm(spacecraft_pos) ** 3
            spacecraft_p = spacecraft_p + F * dt
            spacecraft_pos = spacecraft_pos + (spacecraft_p / spacecraft.mass) * dt
            self.interval.append(spacecraft_pos)
            self.speeds.append(spacecraft_p / spacecraft.mass)
        spacecraft.r = spacecraft_pos
        spacecraft.velocity = spacecraft_p / spacecraft.mass


def goman(spacecraft, target_orbit, start_time):
    if spacecraft.orbit.e != 0:
        print('Не круговые орбиты не поддерживаются')
        return (0, 0)
    if target_orbit.e != 0:
        print('Не круговые орбиты не поддерживаются')
        return (0, 0)
    r1 = spacecraft.orbit.a
    r2 = target_orbit.a
    if r1 < r2:
        k = 1
    else:
        k = -1

    v1 = spacecraft.velocity
    v2 = -sqrt(grav_param / r2) * v1 // np.linalg.norm(v1)
    delta_V_A = sqrt(grav_param / r1) * (sqrt(2 * r2 / (r1 + r2)) - 1) * v1 / np.linalg.norm(v1)
    delta_V_B = sqrt(grav_param / r2) * (1 - sqrt(2 * r1 / (r1 + r2))) * v1 // np.linalg.norm(v1)
    transfer_orbit = Orbit(np.linalg.norm(np.cross(spacecraft.r, v1 + delta_V_A)) ** 2 / grav_param,
                           (r2 - r1) / (r2 + r1))
    # transfer_orbit.draw_orbit()
    T_transfer = pi * sqrt(((0.5 * (r1 + r2)) ** 3) / grav_param)
    spacecraft.velocity = v1 + delta_V_A
    interval_A = Interval(start_time, T_transfer, spacecraft)
    spacecraft.velocity = v2
    spacecraft.orbit = target_orbit
    return (interval_A)

def simulation(track,sc):
    scene.forward = vector(0,0,-1)
    earth = sphere(pos=vector(0,0,0), radius=6.4e6, color=color.blue)
    earth.mass = earth_mass
    colors=[color.red,color.cyan,color.purple,color.green,color.orange,color.white]
    l=0
    for i in track.track_list:
        spacecraft = sphere(pos=vector(i.start_r[0],i.start_r[1],i.start_r[2]), radius=5e5, color=colors[l],
                make_trail=True, interval=10, retain=100000)
        l+=1
        for k in i.interval:
            rate(6000)
            spacecraft.pos = vector(k[0],k[1],k[2])


class track:
    def __init__(self):
        self.track_list = []

    def append_interval(self, interval):
        self.track_list.append(interval)

    def get_track_list(self):
        return (self.track_list)

    def load_from_file(self):
        pass

orbit1=Orbit(20e6,0)
orbit2=Orbit(30e6,0)
orbit3=Orbit(10e6,0)
voyager=Spacecraft(1e4,orbit1)
track1=track()                      #Здесь руками созданы интервалы и итоговый маршрут

interval_1=Interval(0,35000,voyager)
interval_2=goman(voyager,orbit2,20000)
interval_3=Interval(0,35000,voyager)
interval_4=goman(voyager,orbit3,20000)
interval_5=Interval(0,25000,voyager)
voyager.get_momentum(np.array([1e7,0,0]))
interval_6=Interval(0,35000,voyager)
track1.append_interval(interval_1)
track1.append_interval(interval_2)
track1.append_interval(interval_3)
track1.append_interval(interval_4)
track1.append_interval(interval_5)
track1.append_interval(interval_6)

simulation(track1,voyager)          #Запуск функции симуляции



x=[i[0] for i in interval_1.interval]
x1=[i[0] for i in interval_3.interval]
y=[i[1] for i in interval_1.interval]
y1=[i[1] for i in interval_3.interval]
a=[]
q=[]
for i in track1.track_list:
    a.extend(i.speeds)
    q.extend(i.interval)
b=[np.linalg.norm(i) for i in a[::10]]          #Здесь графики на коленке, не очень аккуратные
c_x=[i[0] for i in q]
c_y=[i[1] for i in q]
plot(b)
show()
xlim(-4e7,4e7)
ylim(-4e7,4e7)
plot(c_x,c_y)
show()
#orbit1.draw_orbit()
