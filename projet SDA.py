# -*- coding: utf-8 -*-
"""
Created on Mon Jan 17 07:10:42 2022

@author: ulrich
"""
import numpy as np
import matplotlib.pyplot as plt
class points:
    def __init__(self,masse,x,y,z):
        self.masse, self.x, self.y,self.z = masse,x, y,0
     
    # def __repr__(self):
    #     return '{}: {}'.format(str((self.x, self.y)), repr(self.payload))
    # def __str__(self):
    #     return 'P({:.2f}, {:.2f})'.format(self.x, self.y)
    def distance_to(self, other):
        try:
           other_x, other_y = other.x, other.y
        except AttributeError:
            other_x, other_y = other
        return np.hypot(self.x - other_x, self.y - other_y)
class rectangle:
    def __init__(self, cx, cy, largeur, hauteur,):
        self.cx, self.cy = cx, cy
        self.largeur, self.hauteur = largeur, hauteur
        self.bord_ouest, self.bord_est = cx - largeur/2, cx + largeur/2
        self.bord_nord, self.bord_sud = cy + hauteur/2, cy - hauteur/2
    def __repr__(self):
        return str((self.bord_ouest, self.bord_est,
        self.bord_nord, self.bord_sud))
    def __str__(self):
        return '({:.2f}, {:.2f}, {:.2f}, {:.2f})'.format(self.bord_ouest, self.bord_est,
        self.bord_nord, self.bord_sud)
    def contenir(self, point):
        """Le point (un objet Point ou un tuple (x,y)) est-il à l'intérieur de ce Rect ?"""

        try:
            point_x, point_y = point.x, point.y
        except AttributeError:
            point_x, point_y = point
        return (point_x >= self.bord_ouest and
                point_x <  self.bord_est and
                point_y <= self.bord_nord and
                point_y > self.bord_sud)
    def dessiner(self, ax, c='k', lw=1, **kwargs):
        x1, y1 = self.bord_ouest, self.bord_sud
        x2, y2 = self.bord_est, self.bord_nord
        ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1], c=c, lw=lw, **kwargs)
# =============================================================================
#         plt.axvline(x=self.cx, ymin=self.bord_sud/height, ymax=self.bord_nord/height, c='red')
#         plt.axhline(y=self.cy, xmin=self.bord_ouest/width, xmax=self.bord_est/width, c='red')
# =============================================================================
# =============================================================================
#         ax.plot([x1,x2,x2,x1,x1],[y1,y1,y2,y2,y1], c=c, lw=lw, **kwargs)
# =============================================================================
class quadTree:
    def __init__(self, rectangle, profondeur):
        """ Initialisez ce nœud du quadtree.

        boundary est un objet Rect définissant la région à partir de laquelle les points sont
        placé dans ce nœud ; max_points est le nombre maximum de points que
        le nœud peut tenir avant de devoir se diviser (se ramifier en quatre autres nœuds) ;
        deppth garde une trace de la profondeur de ce nœud dans le quadtree.

        """
        self.rectangle = rectangle
        self.point=None
        self.profondeur = profondeur
        # Un indicateur pour indiquer si ce nœud s'est divisé (ramifié) ou non.
        self.divided = False
        self.nord_ouest=None
        self.sud_ouest=None
        self.nord_est=None
        self.sud_est=None
    def divide(self):
        """Divisez (branchez) ce nœud en générant quatre nœuds enfants."""

     
        # The boundaries of the four children nodes are "northwest",
        # "northeast", "southeast" and "southwest" quadrants within the
        # boundary of the current node.
        self.nord_ouest = quadTree(rectangle((self.rectangle.cx +self.rectangle.largeur/4), (self.rectangle.cy + self.rectangle.hauteur/4), (self.rectangle.largeur/2), (self.rectangle.hauteur/2)),(self.profondeur + 1))
            
        self.sud_ouest = quadTree(rectangle((self.rectangle.cx +self.rectangle.largeur/4), (self.rectangle.cy -self.rectangle.hauteur/4), (self.rectangle.largeur/2), (self.rectangle.hauteur/2)),(self.profondeur + 1))
            
        self.sud_est =quadTree(rectangle((self.rectangle.cx -self.rectangle.largeur/4), (self.rectangle.cy - self.rectangle.hauteur/4), (self.rectangle.largeur/2), (self.rectangle.hauteur/2)),(self.profondeur + 1))
            
        self.nord_est =quadTree(rectangle((self.rectangle.cx -self.rectangle.largeur/4), (self.rectangle.cy + self.rectangle.hauteur/4), (self.rectangle.largeur/2), (self.rectangle.hauteur/2)),(self.profondeur + 1))
            
        self.divided = True
    
    def insert(self, point):
        """Essayez d'insérer un point Point dans ce QuadTree."""

        if not self.rectangle.contenir(point):
            print(self.rectangle.contenir(point))
            print("largeur rectangle=",self.rectangle.largeur," hauteur=", self.rectangle.hauteur)
            print("bord_ouest=",self.rectangle.bord_ouest," bord_est=", self.rectangle.bord_est)
            print("bord_nord=",self.rectangle.bord_nord," bord_sud=", self.rectangle.bord_sud)
            print("point:",point.x," y:",point.y)
            print("le point ne se trouve pas dans l'intervalle")
        
        elif self.divided==False and self.rectangle.contenir(point)==True and self.point is None:
            # s'il n'est pas neccessaire de diviser notre quadtree
            
            
            self.point=point
           
           

        # No room: divide if necessary, then try the sub-quads.
        
        elif not self.divided:
            self.divide()   
            recuperer=points(self.point.masse, self.point.x,self.point.y,self.point.z)
            
            self.point.x=(self.point.masse*self.point.x + point.masse*point.x)/(self.point.masse+point.masse)
            self.point.y=(self.point.masse*self.point.y + point.masse*point.y)/(self.point.masse+point.masse)
            self.point.masse=self.point.masse+point.masse
            self.insert(recuperer)
        else:                      
            cx, cy = self.rectangle.cx, self.rectangle.cy
            largeur, hauteur = self.rectangle.largeur / 2, self.rectangle.hauteur / 2
            if  point.x>self.rectangle.cx and point.y>self.rectangle.cy :
             
                
                self.point.x=(self.point.masse*self.point.x + point.masse*point.x)/(self.point.masse+point.masse)
                self.point.y=(self.point.masse*self.point.y + point.masse*point.y)/(self.point.masse+point.masse)
                self.point.masse=self.point.masse+point.masse
                return self.nord_ouest.insert(point)
            
            elif point.x>=self.rectangle.cx and point.y<=self.rectangle.cy:
               
                self.point.x=(self.point.masse*self.point.x + point.masse*point.x)/(self.point.masse+point.masse)
                self.point.y=(self.point.masse*self.point.y + point.masse*point.y)/(self.point.masse+point.masse)
                self.point.masse=self.point.masse+point.masse
                return self.sud_ouest.insert(point)
            
            elif point.x<self.rectangle.cx and point.y<self.rectangle.cy:
               
                self.point.x=(self.point.masse*self.point.x + point.masse*point.x)/(self.point.masse+point.masse)
                self.point.y=(self.point.masse*self.point.y + point.masse*point.y)/(self.point.masse+point.masse)
                self.point.masse=self.point.masse+point.masse
                return self.sud_est.insert(point)
            else:
                
                self.point.x=(self.point.masse*self.point.x + point.masse*point.x)/(self.point.masse+point.masse)
                self.point.y=(self.point.masse*self.point.y + point.masse*point.y)/(self.point.masse+point.masse)
                self.point.masse=self.point.masse+point.masse
                return self.nord_est.insert(point)    

    

    def rechercher_point(self,point):
        if not self.rectangle.contenir(point):
            # If the domain of this node does not intersect the search
            # region, we don't need to look in it for points.
            return False

        if not self.divided and self.point!=None:
            if point.x==self.point.x and point.y==self.point.y:
                return True
            else:
                return False
        if self.divided:                  
            if point.x<self.rectangle.bord_est and point.x>self.rectangle.cx and point.y>self.rectangle.cy and point.y<self.rectangle.bord_nord:
                print("nord_ouest")
                return self.nord_ouest.rechercher_point(point)
            
            elif point.x<=self.rectangle.bord_est and point.x>=self.rectangle.cx and point.y<=self.rectangle.cy and point.y>=self.rectangle.bord_sud:
                print("sud_ouest")
                return self.sud_ouest.rechercher_point(point)
            
            elif point.x<self.rectangle.cx and point.x>self.rectangle.bord_ouest and point.y<self.rectangle.cy and point.y>self.rectangle.bord_sud:
                print("sud_est")
                return self.sud_est.rechercher_point(point)
            else:
                print("nord_est")
                return self.nord_est.rechercher_point(point)
# =============================================================================
#     def rechercher_point_delimiter(self,point_trouver,point):
#          if not self.rectangle.contenir(point):
#             # If the domain of this node does not intersect the search
#             # region, we don't need to look in it for points.
#             return False
# 
#          if not self.divided :
#             if point.x==self.point.x and point.y==self.point.y:
#                 point_trouver.append(point)
#             else:
#                 return False
#          
#       
#          if point.x<self.rectangle.bord_est and point.x>self.rectangle.cx and point.y>self.rectangle.cy and point.y<self.rectangle.bord_nord:
#             print("nord_ouest")
#             return self.nord_ouest.rechercher_point(point)
#         
#          elif point.x<=self.rectangle.bord_est and point.x>=self.rectangle.cx and point.y<=self.rectangle.cy and point.y>=self.rectangle.bord_sud:
#             print("sud_ouest")
#             return self.sud_ouest.rechercher_point(point)
#         
#          elif point.x<self.rectangle.cx and point.x>self.rectangle.bord_ouest and point.y<self.rectangle.cy and point.y>self.rectangle.bord_sud:
#             print("sud_est")
#             return self.sud_est.rechercher_point(point)
#          else:
#             print("nord_est")
#             return self.nord_est.rechercher_point(point)
# =============================================================================


    def affichernoeudExterne(self):
        if self.divided==False and self.point is not None:
            print("point X=",self.point.x," Y=",self.point.y, "profondeur=",self.profondeur, "masse=",self.point.masse)
        else:
            
            if self.nord_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.nord_ouest.affichernoeudExterne()
            if self.sud_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.sud_ouest.affichernoeudExterne()
            if self.nord_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.nord_est.affichernoeudExterne()
            if self.sud_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.sud_est.affichernoeudExterne()
    def afficherNoeurdInterne(self):
        if self.divided==True and self.point is not None:
            print("point X=",self.point.x," Y=",self.point.y, "profondeur=",self.profondeur, "masse=",self.point.masse)
            if self.nord_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.nord_ouest.afficherNoeurdInterne()
            if self.sud_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.sud_ouest.afficherNoeurdInterne()
            if self.nord_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.nord_est.afficherNoeurdInterne()
            if self.sud_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.sud_est.afficherNoeurdInterne()
    def draw(self, ax):
        """Draw a representation of the quadtree on Matplotlib Axes ax."""
        if self.divided==True:
            print("point: X=",self.point.x," y=",self.point.y)
        self.rectangle.dessiner(ax)
        if self.divided:
            self.nord_ouest.draw(ax)
            self.nord_est.draw(ax)
            self.sud_est.draw(ax)
            self.sud_ouest.draw(ax)
    def recuperer_point(self,liste):
        if self.divided==False and self.point is not None:
            liste.append(self.point)
     
        if self.sud_est!=None:
               print("toucher_sud_est")
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.sud_est.recuperer_point(liste)
        
               
        if self.sud_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               print("toucher surd_ouest")
               self.sud_ouest.recuperer_point(liste)
        if self.nord_est!=None:
               print("toucher nord_est ")
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               self.nord_est.recuperer_point(liste)
        if self.nord_ouest!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
               print("toucher nord_ouest")
               self.nord_ouest.recuperer_point(liste)

    def calculer_force(self,quadtree,point,liste,listeForce,g):
        if quadtree.divided==False and quadtree.point is not None:
            if point.x!=quadtree.point.x and point.y!=quadtree.point.y:
                F=(g*point.masse*quadtree.point.masse)/(point.distance_to(quadtree.point)**2)
                liste.append(points(quadtree.point.masse,quadtree.point.x,quadtree.point.y,quadtree.point.z))
                listeForce.append(F)
        elif quadtree.divided==True and quadtree.point is not None:
            d=point.distance_to(quadtree.point)
            if d>2*quadtree.rectangle.hauteur:
                F=(g*point.masse*quadtree.point.masse)/(point.distance_to(quadtree.point)**2)
                liste.append(points(quadtree.point.masse,quadtree.point.x,quadtree.point.y,quadtree.point.z))
                listeForce.append(F)
            else:
                if quadtree.nord_ouest!=None:
                    self.calculer_force(quadtree.nord_ouest, point,liste,listeForce, g)
                if quadtree.sud_ouest!=None:
                    self.calculer_force(quadtree.sud_ouest, point,liste,listeForce, g)
                if quadtree.nord_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
                    self.calculer_force(quadtree.nord_est, point,liste,listeForce, g)
                if quadtree.sud_est!=None:
               ##print("rectangle: masse=",self.point.masse," X=",self.point.x," Y=",self.point.y)
                    self.calculer_force(quadtree.sud_est, point,liste, listeForce,g)
                
            

            
        
        
        
        
        
        
    


def getAcc( pos, mass, G, softening ):
	"""
    Calculer l'accélération de chaque particule grâce à la loi de Newton
pos est une matrice N x 3 de positions
la masse est un vecteur N x 1 de masses
G est la constante gravitationnelle de Newton
ramollissement est la longueur de ramollissement
a est la matrice N x 3 des accélérations
	"""
    
	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r^3 for all particle pairwise particle separations 
	inv_r3 = (dx**2 + dy**2 + dz**2 + softening**2)
	inv_r3[inv_r3>0] = inv_r3[inv_r3>0]**(-1.5)

	ax = G * (dx * inv_r3) @ mass
	ay = G * (dy * inv_r3) @ mass
	az = G * (dz * inv_r3) @ mass
	
	# pack together the acceleration components
	a = np.hstack((ax,ay,az))

	return a   
def getEnergy( pos, vel, mass, G ):
	"""
	Get kinetic energy (KE) and potential energy (PE) of simulation
	pos is N x 3 matrix of positions
	vel is N x 3 matrix of velocities
	mass is an N x 1 vector of masses
	G is Newton's Gravitational constant
	KE is the kinetic energy of the system
	PE is the potential energy of the system
	"""
	# Kinetic Energy:
	KE = 0.5 * np.sum(np.sum( mass * vel**2 ))


	# Potential Energy:

	# positions r = [x,y,z] for all particles
	x = pos[:,0:1]
	y = pos[:,1:2]
	z = pos[:,2:3]

	# matrix that stores all pairwise particle separations: r_j - r_i
	dx = x.T - x
	dy = y.T - y
	dz = z.T - z

	# matrix that stores 1/r for all particle pairwise particle separations 
	inv_r = np.sqrt(dx**2 + dy**2 + dz**2)
	inv_r[inv_r>0] = 1.0/inv_r[inv_r>0]

	# sum over upper triangle, to count each interaction only once
	PE = G * np.sum(np.sum(np.triu(-(mass*mass.T)*inv_r,1)))
	
	return KE, PE;  




 	
 	
# =============================================================================
#  	# Simulation parameters
# N         = 50    # Number of particles
# t         = 0      # current time of the simulation
# tEnd      = 10.0   # time at which simulation ends
# dt        = 0.01   # timestep
# softening = 0.1    # softening length
# G         = 1.0    # Newton's Gravitational Constant
# plotRealTime = True # switch on for plotting as the simulation goes along
#  	
#  	# Generate Initial Conditions
# np.random.seed(17)            # set the random number generator seed
# quatree=quadTree(rectangle(5, 5, 10, 10),0)
#     
#coordx = np.random.randn(N, 2) * height/3 + (width/2, height/2)
#coordy= np.random.randn(N, 2) * height/3 + (width/2, height/2)
#i=0
#coords=[]
#print((np.shape(coordx))[0])

#while i<N:
#    point=points(coordx[i][0], coordx[i][1], coordy[i][0], coordy[i][1])
#    print(point)
#    coords.append(point)
#    i=i+1

# mass = coords[:][0]
#     # total mass of particles is 20
# pos  = np.random.randn(N,3)
#pos[:][0]=coords[:][1]
#pos[:][1]=coords[:][2]
#pos[:][2]=coords[:][3]
# i=0
# while (i< N):
#         X=pos[:][0][:]
#         
#         point=points(mass[i],pos[i][0],pos[i][1],pos[i][2])
#         quatree.insert(point)
# 
#         # randomly selected positions and velocities
# vel  = np.random.randn(N,3)
# listepoint=[]
# quatree.recuperer_point(listepoint)
# # Convert to Center-of-Mass frame
# vel -= np.mean(mass * vel,0) / np.mean(mass)
#  	
#  	# calculate initial gravitational accelerations
# acc = getAcc( pos, mass, G, softening )
#  	
#  	# calculate initial energy of system
# KE, PE  = getEnergy( pos, vel, mass, G )
#  	
#  	# number of timesteps
# Nt = int(np.ceil(tEnd/dt))
#  	
#  	# save energies, particle orbits for plotting trails
# pos_save = np.zeros((N,3,Nt+1))
# pos_save[:,:,0] = pos
# KE_save = np.zeros(Nt+1)
# KE_save[0] = KE
# PE_save = np.zeros(Nt+1)
# PE_save[0] = PE
# t_all = np.arange(Nt+1)*dt
#  	
#  	# prep figure
# fig = plt.figure(figsize=(4,5), dpi=80)
# grid = plt.GridSpec(3, 1, wspace=0.0, hspace=0.3)
# ax1 = plt.subplot(grid[0:2,0])
# ax2 = plt.subplot(grid[2,0])
#  	
#  	# Simulation Main Loop
# for i in range(Nt):
#     
# 		# (1/2) kick
# 	vel += acc * dt/2.0
# 		
# 		# drift
# 	pos += vel * dt
# 		
# 		# update accelerations
# 	acc = getAcc( pos, mass, G, softening )
# 		
# 		# (1/2) kick
# 	vel += acc * dt/2.0
# 		
# 		# update time
# 	t += dt
# 		
# 		# get energy of system
# 	KE, PE  = getEnergy( pos, vel, mass, G )
# 		
# 		# save energies, positions for plotting trail
# 	pos_save[:,:,i+1] = pos
# 	KE_save[i+1] = KE
# 	PE_save[i+1] = PE
# 		
# 		# plot in real time
# if plotRealTime or (i == Nt-1):
#  	plt.sca(ax1)
#  	plt.cla()
#  	xx = pos_save[:,0,max(i-50,0):i+1]
#  	yy = pos_save[:,1,max(i-50,0):i+1]
#  	plt.scatter(xx,yy,s=1,color=[.7,.7,1])
#  	plt.scatter(pos[:,0],pos[:,1],s=10,color='blue')
#  	ax1.set(xlim=(-2, 2), ylim=(-2, 2))
#  	ax1.set_aspect('equal', 'box')
#  	ax1.set_xticks([-2,-1,0,1,2])
#  	ax1.set_yticks([-2,-1,0,1,2])
#  			
#  	plt.sca(ax2)
#  	plt.cla()
#  	plt.scatter(t_all,KE_save,color='red',s=1,label='KE' if i == Nt-1 else "")
#  	plt.scatter(t_all,PE_save,color='blue',s=1,label='PE' if i == Nt-1 else "")
#  	plt.scatter(t_all,KE_save+PE_save,color='black',s=1,label='Etot' if i == Nt-1 else "")
#  	ax2.set(xlim=(0, tEnd), ylim=(-300, 300))
#  	ax2.set_aspect(0.007)
#  			
#  	plt.pause(0.001)
#  	    
#  	
#  	
#  	# add labels/legend
# plt.sca(ax2)
# plt.xlabel('time')
# plt.ylabel('energy')
# ax2.legend(loc='upper right')
#  	
#  	# Save figure
# plt.savefig('nbody.png',dpi=240)
# plt.show()
#  	    
# 
# =============================================================================
 	


  
   




DPI = 72
np.random.seed(60)

width, height = 600, 400

N = 1000000
coordx = np.random.randn(N, 2) * height/3 + (width/2, height/2)
coordy= np.random.randn(N, 2) * height/3 + (width/2, height/2)
i=0
coords=[]
print((np.shape(coordx))[0])

while i<N:
    point=points(coordx[i][0], coordx[i][1], coordy[i][0], coordy[i][1])
    print(point)
    coords.append(point)
    i=i+1

print("-----",coords)



domain = rectangle(width/2, height/2, width, height)
qtree = quadTree(domain, 0)

for point in coords:
    qtree.insert(point)
# =============================================================================
# print(qtree.__str__())
# print('Nombre de points dans le domaine =', len(qtree))
# =============================================================================

fig = plt.figure(figsize=(700/DPI, 500/DPI), dpi=DPI)
ax = plt.subplot()
ax.set_xlim(0, width)
ax.set_ylim(0, height)
qtree.draw(ax)

ax.scatter([p.x for p in coords], [p.y for p in coords], s=4)

ax.set_xticks([])
ax.set_yticks([])

# =============================================================================
# centre, rayon = (width/2, height/2), 120
# found_points = []
# qtree.query_radius(centre, rayon, found_points)
# print('Number of found points =', len(found_points))
# 
# ax.scatter([p.x for p in found_points], [p.y for p in found_points],
#            facecolors='none', edgecolors='r', s=32)
# 
# circle = plt.Circle(centre, rayon, ec='r')
# rectangle(*centre, 2*rayon, 2*rayon).draw(ax, c='r')
# =============================================================================

ax.invert_yaxis()
plt.tight_layout()
plt.savefig('search-quadtree-circle.png')
plt.show()





rect=rectangle(5,5,10,10)          
application=quadTree(rect, 0)
p=points(2, 2, 2,2)
p1=points(4,2,7,7)
p2=points(10,7,7,7)
p3=points(10,7,2,2)
P4=points(4,4,4,4)
P5=points(1,3,1,2)
P6=points(1,3,0.5,2)
print(rect.contenir(p))

application.insert( point=p)
application.insert(point=p1)
application.insert(point=p2)
application.insert(point=p3)
application.insert(point=P4)
application.insert(point=P6)
application.insert(point=P5)
liste=[]
liste_point=application.recuperer_point(liste)
for i in liste:
    print("liste:",i)

print("recherche de P5", application.rechercher_point(p3))
print("----------noeud externe--------------------")
application.affichernoeudExterne()
force=[]
listenoeud=[]

application.calculer_force(application, p, listenoeud,force, 1)
for i in listenoeud:
    print("noeud:" ,i.x,"----",i.y,i.masse)
for j in force:
    print("force:",j)

print(application.afficherNoeurdInterne())




