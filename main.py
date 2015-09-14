import numpy as np
import math
from collections import namedtuple

# Tables: http://www.openelectrical.org/wiki/index.php?title=AACSR_Conductors

conductorSpacing = namedtuple("conductorSpacing", "Dab Dbc Dca")
conductorStruct = namedtuple("conductorStruct", "constructionType metalType numCoresTot numSegments driedAndImpreg insulationType Spacing")

class sequence_impedances_cables():
    def __init__(self, f=50, conductor):
        self.w = 2 * np.pi * f # system frequency
        self.u_0 = 4 * np.pi * 10 ** (-7) # permittivity of free space (constant)
        self.conductor = conductor # conductor descriptor
        self.S # Conductor's cross-sectional area (mm^2)
        self.d_c = 2 * (self.S / np.pi) ** 0.5 # Conductor diameter (mm)
        self.R20 # Conductors resistance at 20 deg. C
        self.a20 # Conductors thermal coeffcient at 20 deg. C
        self.t # Actual temperature of the conductor
        
        self.determineKconst(conductor)
        
    def determineKconst(self, conductor):
        defaultList = ["circular", "stranded", "compacted", "sectored"]
        specialList = ["segmental", "segment"]
        if conductor.constructionType in defaultList:
            self.k_s = 1
            self.k_p = 1
        if conductor.driedimpreg and conductor.metal is "Copper":
            self.k_p = 0.8
        if conductor.constructiontype in specialList:
            self.k_p = 0.37
            if conductor.metal is "Copper":
                self.k_s = 0.435
            if conductor.metal is "Aluminium":
                if conductor.segmentnum == 4:
                    self.k_s = 0.28
                elif conductor.segmentnum == 5:
                    self.k_s = 0.19
                elif conductor.segmentnum == 6:
                    self.k_s = 0.12
        
    def _reactance_L(self, w, L): 
        '''Find X_L'''
        return w * L
    
    def _reactance_C(self, w, C):
        '''Find X_C'''
        return 1 / (w * C)
    
    def dc_resistance(self):
        # IEC 60287-1 Clause 2.1
        '''S: cross sectional area in mm^2
        R20: conductor resistance at 20 deg. C
        a20: temperature coeffceient of the material at 20 deg. C
        t: temperature of conductor (deg. C)
        Returns the resistance at temperature t'''
        # DC current assumes no skin and proximity effect, so the conductor is a uniform temperature
        return self.R20 * 1.02 * (10 ** 6) / self.S * (1 + self.a20 * (self.t - 20)) # The conductors resistance in ohms at temp. t in deg. C
    
    def ac_resistance(self, conductor):
        # IEC 60287-1 Clause 2.1
        '''y_s: skin effect factor
        y_p: proximity effect factor
        '''
        y_s = self._skin_factor()
        y_p = self._proximity_factor(conductor)
        R_ac = self.dc_resistance() * (1 + y_s + y_p)
        return R_ac
    
    def _skin_factor(self):
        '''Determine the skin effect factor'''
        dc_res = self.dc_resistance()
        x_4_s = self.x_4_s(self.k_s, dc_res)
        y_s = x_4_s / (192 + x_4_s) #y_s = x_4_s / (192 + 0.8 * x_4_s)
        return y_s
    
    def _proximity_factor(self, conductor):
        '''Determine the proximity effect factor'''
        s = 1 # TODO: change this
        dc_res = self.dc_resistance()
        x_4_s = self.x_4_s(self.k_p, dc_res)
        factor1 = x_4_s / (192 + 0.8 * x_4_s) * (self.d_c/s) ** 2 # s is the distance between conductor axes (mm)
        
        if conductor.numCoresTot == 2:
            factor2 = 2.9
        elif conductor.numCoresTot == 3:
            factor2 = 0.312 * (self.d_c/s) ** 2 + 1.18 / (0.27 + x_4_s / (192 + 0.8 * x_4_s))
        return factor1 * factor2
    #http://www.openelectrical.org/wiki/index.php?title=Cable_Impedance_Calculations 
    
    def x_4_s(self, k, dc_res):
        '''k: constant Ks or Kp'''
        x_4_s = (8 * np.pi * self.f / dc_res * k * 10 ** (-7)) ** 2
        assert x_4_s ** 1/4.0 <= 2.8, "Xs is beyond the limit for an accurate calculation"
        return x_4_s
        

# The self-impedance requires GMR
# The mutual-impedance requires GMD
class overheadLines():
    # https://psspy.org/psse-help-forum/question/94/how-do-you-determine-the-line-constants-for-entering-the-branch-details/
    def __init__(self, f):
        # Constants
        self.w = 2.0 * np.pi * f # system frequency
        self.u_0 = 4.0 * np.pi * 10 ** (-7.0) # permittivity of free space
        self.c = 3.0e8 # speed of light
        self.e_0 = 1.0 / (self.u_0 * self.c**2) # Epsilon 0
        self.p = 100.0 # Earth resistivity in ohmmeters/169.84 ohmmeters 
        self.D_e = 658.37 * (self.p / self.f) ** 0.5 # Eqivilant depth
        
        # Conductor Spacing
        self.Dab
        self.Dbc
        self.Dca
    
        def pos_seq(self, Dab, Dbc, Dca, r_c):
            '''Dab - Dbc - Dca are the spacing between the conductors in meters. 
            GMRc is the geometric mean radius of the conductor'''
            '''Positive (and negative) sequence impedance of a three phase conductor
            '''
            GMD = self._GMD(Dab, Dbc, Dca)
            
            # Assuming ia + ib + ic = 0
            L_a = 2 * 10**(-7) * np.log(GMD/self.GMRc) # L_a = 2 * 10**(-7) * np.log(GMD/GMRu) henry/meter?     
            X_L = self._reactance_L(self.w, L_a)
            Z_1 = r_c + 1j * (X_L) # r_c is the conductors real resistance
            
            #Z_1  = r_c + j * 4 * np.pi * f * 10**(-7) * np.log(GMD/GMRc)
            return Z_1 # ohms/meter
        
        def zero_seq(r_p, GMRc):
            '''r_p: real resistance
            returns the ohms/m'''
            #Z_self = r_p + 9.869e-1 * self.f + 1j * (4e-1 * np.pi * self.f * np.log(self.D_e / GMRc)) # GMRc is the geometric mean radius of a single conductor
            #Z_mutual = 9.869e-1 * self.f + 1j * (4e-1 * np.pi * self.f * np.log(self.D_e / Dpq)) # Dpq is the distance between conductors p and q
            Z_u = r_p / 3 + 9.869e-1 * self.f + 1j * (4e-1 * np.pi * self.f * np.log(self.D_e / self._GMRu(GMRc)))
            Z_0 = 3 * Z_u
            return Z_0
            
        
        def _GMD1(self, *args):
            # Follows the EEA book
            '''Find the equivalent distance also known as the
            Geometric Mean Distance of the phase conductor system
            (assumes no earth wires - will cause an error)'''
            
            assert len(args) % 3 == 0, "Not a valid 3 phase system. Error 1: Integer multiple."
            n = (np.log(len(args) / 3) / np.log(3))
            assert n.is_integer(), "Not a valid 3 phase system. Error 2: Geometric sequence."
            assert len(args) == 3**n, "Not a valid 3 phase system. Error 3: Incorrect number of distances supplied."
            
            radicandTot = 1
            for radicand in args:
                radicandTot *= radicand
            GMD = (radicandTot) ** (1.0 / len(args))
            return GMD
        
        def _GMD2(self, *args):
            # Follows the skm-eleksys calculations, the args will be diffrent to GMR1()
            '''Find the equivalent distance also known as the
            Geometric Mean Distance of the phase conductor system
            (assumes no earth wires - will cause an error)'''
            
            #assert len(args) % 3 == 0, "Not a valid 3 phase system. Error 1: Integer multiple."
            n = (np.log(len(args) / 3) / np.log(3))
            #assert n.is_integer(), "Not a valid 3 phase system. Error 2: Geometric sequence."
            #assert len(args) == 3**n, "Not a valid 3 phase system. Error 3: Incorrect number of distances supplied."
            
            radicandTot = 1
            for radicand in args:
                radicandTot *= radicand
            GMD = (radicandTot) ** (1.0 / len(args))
            return GMD
        
        
        def _GMRu(self, GMRc, num_circuits=1, *args):
            '''Ficticious conductor equivilant of 3 sparse 3-phase conductors
            Geometric Mean Radius
            n_c: is the number of bundled conductors'''
            index = (3.0 * num_circuits) ** 2 # 3.0 means a 3 phase system i.e. 3 conductors
            radicandTot = GMRc ** (3.0 * num_circuits)
            for radicand in args:
                radicandTot *= radicand
            GMRu = radicandTot ** (1.0 / index) # Radical/root
            #GMRu = (GMRc**3 * self.Dab**2 * self.Dbc**2 * self.Dca**2) ** (1 / 9.0) # for a single circuit, 3-phase, 1 conductor per phase
            
            # NB: the real resistance of this ficticious conductor is r_c/3
            return GMRu

class carson_equations(overheadLines):
    '''Carson equations only deal with zero-sequence self and mutual impdances'''
    def __init__(self):
        pass
    
    def self_Z0(self):
        # Assuming this is in ohms/km
        Z11 = rii + (10e-4 * 8 * np.pi * self.f * Pii) + 1j * (8 * np.pi * self.f * 10e-4) * (0.5 * np.log(2 * hi / GMRii) + Qii)
    
    def mutual_Z0(self):
        Z12 = (10e-4 * 8 * np.pi * self.f * Pij) + 1j * (8 * np.pi * self.f * 10e-4) * (0.5 * np.log(Sij / sij) + Qij)


class maxwell_cap(overheadLines):
    '''A class to calculate the parrallel capacitance of an overhead line'''
    def __init__(self):
        pass
    
    def C_matrix(self, P):
        '''P: matrix of potential coeffcients (nxn), 
        n = number of conductors in the system'''
        #C = P ** (-1) # requires a sympy matrix object
        C = np.linalg.inv(P)
        return C
    
    def P_matrix(self):
        const = 1.0 / (2.0 * np.pi * self.e_0)


# Class structures (data containers)
class constants():
    def __init__(self, f):
        '''Electrical constants. Some are frequency dependent.'''
        self.f = f                              # System frequency (Hz)
        self.w = 2.0 * np.pi * self.f                # System frequency (rad/s)
        self.u_0 = 4.0 * np.pi * 10 ** (-7.0)   # Mu 0 (Vacuum permeability) (N/A^2)
        self.c = 3.0e8                          # Speed of light (m/s)
        self.e_0 = 1.0 / (self.u_0 * self.c**2) # Epsilon 0 (permittivity of free space) (F/m)
        self.p = 100.0 #169.84                  # Earth resistivity (ohm meters)
        self.D_e = 658.37 * (self.p / self.f) ** 0.5 # Equivalent  depth (m)


class conductors(constants):
    def __init__(self, n, GMRc=None):
        '''A container class for the physical properties of the conductor itself.
        n: number of conductors comprising this bundle (typically: 1, 2, 3, or 4)'''
        # Get GMRc/D_s/self-GMD form a datasheet preferably (ACSR conductors)
        if GMRc is None:
            self.GMR_c = self._GMRc(r, u_wire, n, R)
        else:    
            self.GMR_c = GMRc
    
    # The GMR of the conductor itself
    def _GMRc(self, r, u_wire, n=1, R=0):
            '''
            n: number of conductor strands in a single phase conductor
            r: radius of a individual strand of conductor
            R: outside radius of stranded conductor
            Although the formulas in this module are derived considering transposition, the 
            same formulas are also used for non-transposed cases to get approximate values.
            '''
            # https://uqu.edu.sa/files2/tiny_mce/plugins/filemanager/files/4300303/Transmission%20of%20EE/t2%2009%2010/Transmission%20and%20Distribution%20of%20Electrical%20Power%20-%203%20Electric%20Transmission%20Line%20Parameters.pdf
            # https://en.wikipedia.org/wiki/Permeability_(electromagnetism)
            # http://physics.stackexchange.com/questions/10827/is-aluminium-magnetic
            u_r = u_wire / self.u_0
            alpha = 1.0 / 4.0 * u_r # alpha ~= 1/4 if a material is non-magnetic or paramagnetic (e.g. Al)
            r1 = r * np.exp(-alpha)
            radicand = 1.0
            
            if n == 1:
                radicand = r1 * (2.0*r)**(n-1)
            elif n == 2:
                radicand = r1 * (2.0*r)**(n-1)
            elif n == 3:
                radicand = r1 * (2.0*r)**(n-1)
            elif n == 4:
                radicand = r1 * (2.0*r)**(n-1) * (2)**(1.0/2.0)
            elif n == 7:
                radicand = r1 * (2.0*r)**(n-1) * (2)**(6.0/7.0) * (3)**(6.0/7.0)
                # should be the same as: GMRc = 0.726 * R

            # Determine the irrational number from the radical expression
            index = n
            GMRc = radicand ** (1.0 / index)
            
            # I can not derive the general equation, so these values are hard-coded
            # https://books.google.co.nz/books?id=xgSXRWvVwaEC&pg=PA12&lpg=PA12&dq=self-GMD+or+GMR+of+stranded+conductors&source=bl&ots=90g4DEVpeJ&sig=NSGVU_7X9k4lCM_mq5m1QbD5358&hl=en&sa=X&ved=0CB0Q6AEwAGoVChMI6_GH1ovcxwIVAxmmCh1o3QDx#v=onepage&q=self-GMD%20or%20GMR%20of%20stranded%20conductors&f=false
            if n == 19:
                GMRc = 0.758 * R
            elif n == 37:
                GMRc = 0.768 * R
            elif n == 61:
                GMRc = 0.772 * R
            elif n == 91:
                GMRc = 0.774 * R
            elif n == 127:
                GMRc = 0.776 * R
                
            # Hollow standed and ACSR conductors (neglects steel stands)
            if n == 30:
                GMRc = 0.826 * R
            elif n == 26:
                GMRc = 0.809 * R
            elif n == 54:
                GMRc = 0.810 * R
            
            return GMRc
        

class circuit(conductors, constants): 
    def __init__(self, earthWire=(0,0), nb=(1, 0), *args):
        '''A container class/struct(ure) for the electrical circuit.
        earthWire: the location of an earth wire (if present), leave (0,0) if there are no earth wires
        nb: tuple that represent (the number of conductors in each phase, spacing between conductors)
        *args: n number of tuple arguments that represent the coordinates (x, y) of the phase conductors 
        '''
        '''
        Pass in n tuples of coordinates (x, y) that represent the spacing of the phase conductors - max. 1 circuit input.
        All dimensions are relative the bottom center of a pole (2D plane only/cross-section).
        '''
        if earthWire == (0,0):
            self.eWire = None # earth wire
        else:
            self.eWire = earthWire
        self.nb = nb
        self.coords = [xy for xy in args] # unpack the (x, y) args
    
    # Calculate GMDs
    
    # GMR of multi conductor per phase system (aka "bundled conductors")
    # Note that the use of self.coords is not ideal, as no-one will enter coords for a each condcutor in one phase
    def GMR(self, n, r, d=1):
        '''
        n: number of conductors in the bundle
        r: equal to GMR_c = D_s
        d: distance between conductors of the same phase (that form a bundle). Assign 1 for a single conductor bundle.
        
        This function returns the GMR for a bundled group of conductors (i.e. 1 phase)
        The function dimensions from every conductor (node) to every other conductor (node),
        such that no single distance between conductors is dimensioned more than once. 
        
        This function is only valid for bundles up to 4 conductors large - needs coorected coord spacing.
        '''
        nodeHistory = []
        distances = []
        # Check every node in the bundle
        for node in self.coords:
            nodeHistory.append(node)
            (x1 , y1) = node
            # For this particular node, find all other nodes in the bundle
            # that are not this node itself, or haven't already been dimentioned
            # releative to this node.
            for node1 in self.coords:
                if node1 not in nodeHistory:
                    (x2, y2) = node1
                    distances.append(math.hypot(x2-x1, y2-y1))
        assert len(distances) == n*(n-1)/2, "GMR error: Not a triangular number"
        # This has been simplified by assuming all conductors and distances (d) are identical
        exponent = 1.0 / len(distances) # In full form: n * (1/n^2) = 1/n
        radicand = r # This ensures that even for a single conductor the right GMR is returned
        for d in distances:
            radicand *= d
        GMR = radicand ** exponent
        return GMR
    
    def GMR_gen_coords(self, n, d=0):
        #http://math.stackexchange.com/questions/1267646/finding-vertices-of-regular-polygon
        '''
        Generate as set of coordinates for n bundled conductors
        seperateed by distance d (i.e. length of a 1 side of a polygon)
        '''
        pt1 = (0, 0)
        pt2 = (d, 0)
        coords = []
        #assert n >= 1, "Not a valid number of conductors"
        #assert isinstance(n, int), "Not a valid number of conductors (float)"
        #assert d >= 0, "Not a valid distance between conductors"
        #assert n >= 2 and d <= 0, "More than 1 conductor. Specify a seperation distance" 
        # Single conductor case
        if n == 1 or d <= 0:
            coords += [pt1]
        # Dual conductor case
        elif n == 2:
            coords += [pt1, pt2]
        # n sided polygon case (>2 conductors)
        elif n >= 3:
            v0 = (d/2.0, d/( 2.0 * np.tan(np.pi/n))) # The centre point of the polygon (between pt1 and pt2)
            coords = coords + [tuple(np.subtract(pt1, v0))] + [tuple(np.subtract(pt2, v0))] # Offset the known edge/side coords
            theta = 2.0 * np.pi / n # The angle two adjacent verices make at the origin
            R_theta = [[np.cos(theta), -np.sin(theta)], [np.sin(theta), np.cos(theta)]] # The 2D rotation matrix
            for side in range(n-2): # Minus 2 because we already have two points
                vertex = np.dot(R_theta, coords[side+1])
                coords.append(tuple(vertex))
            # This returns the coords to the original origin point - can be removed since we only care about GMD/GMR
            for c in range(len(coords)):
                coords[c] = tuple(np.add(v0, coords[c]))
        return coords


class arrangement(circuit, constants):
    def __init__(self, conductor, *args):
        # remove condcutor from the __init__ this should be accessable from class inheritncae
        '''
        conductor: is class/struct/container of the all the physical properties of
        the conductor. 
        args: represents n number of "circuit" objects that will combined to form an overall arrangement
        composed of multiple circuits, earth wires.'''
        self.conductor = conductor
        # combine n circuits to form a total arrangement of conductors
        self.conductorGroup = []
        for circuit in args:
            self.conductorGroup.append(circuit.coords)
    
    # The GMD of a poly-phase group of conductors (aka a circuit/3 phase system)
    def GMD(self, a):
        distances = []
        for i in range(len(phaseGroup1)):
            (x1 , y1) = phaseGroup1[i]
            for j in range(len(phaseGroup2)):
                if i != j:
                    (x2, y2) = phaseGroup2[j]
                    distances.append(((x1-x2)**2 + (y1-y2)**2)**0.5)
        assert len(distances) != 0, "Not a complete circuit. No support for SWER systems."
        radicad = 1
        for d in distances:
            radicand *= d
        GMD = radicand ** (1.0/len(distances))
        return GMD
    
    def GMD(self, n, r, d=1):
        '''
        n: number of (equivalent) phase conductors
        r: equal to GMR_c = D_s
        d: distance between (equivalent) conductors of a different phase
        
        This function returns the GMD for a bundled group of conductors (i.e. 1 phase)
        The function dimensions from every conductor (node) to every other conductor (node),
        such that no single distance between conductors is dimensioned more than once. 
        
        This function is only valid for bundles up to 4 conductors large.
        '''
        nodeHistory = []
        distances = []
        # Check every conductor (node) in the poly-phase system
        for node in self.coords:
            nodeHistory.append(node)
            (x1 , y1) = node
            # For this particular conductor (node), find distances to all phase conductors (nodes)
            # that are not this conductor (node) itself, its transpose conductor (node) (double circuit only)
            # or any distance path(s) already dimentioned.
            for node1 in self.coords:
                if node1 not in nodeHistory:
                    (x2, y2) = node1
                    distances.append(math.hypot(x2-x1, y2-y1))
        assert len(distances) == n*(n-1)/2, "GMR error: Not a triangular number"
        # This has been simplified by assuming all conductors and distances (d) are identical
        exponent = 1.0 / len(distances) # In full form: n * (1/n^2) = 1/n
        radicand = r # This ensures that even for a single conductor the right GMR is returned
        for d in distances:
            radicand *= d
        GMR = radicand ** exponent
        return GMR
    
from sympy.solvers import solve
from sympy import symbol

def carson_infinite_seriesA(mode, *coords):
    # unpack *coords
    
    if mode == "self":
        theta = 0 # sin^-1(x_ij)/S_ij
        k = 10e-3 * 4
    elif mode == "mutual":
        pass
    else:
        assert 0, "Not valid input"
    return theta, k

def carson_ininite_seriesB(theta, k):
    '''The simplfied implmentation'''
    P = np.pi / 8
    Q = -0.0386 + 0.5*np.log(2/k)
    return (P, Q)
