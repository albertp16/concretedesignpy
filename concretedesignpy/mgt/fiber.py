import math

def area_diam(diameter: float) -> float:
    """
    Return the cross-sectional area of a circle given its diameter.

    Args:
        diameter (float): The diameter of the circle (must be positive).

    Returns:
        float: The area of the circle.

    Raises:
        ValueError: If 'diameter' is not positive.
    """
    if diameter <= 0:
        raise ValueError("Diameter must be positive.")
    return (math.pi / 4) * (diameter ** 2)


def steel_area(num_bars: int, bar_area: float) -> float:
    """
    Return total steel area for a given number of bars.

    Args:
        num_bars (int): Number of bars (must be positive).
        bar_area (float): Cross-sectional area of one bar (must be positive).

    Returns:
        float: The total steel area.

    Raises:
        ValueError: If 'num_bars' <= 0 or 'bar_area' <= 0.
    """
    if num_bars <= 0:
        raise ValueError("Number of bars must be positive.")
    if bar_area <= 0:
        raise ValueError("Bar area must be positive.")
    return num_bars * bar_area


def area_ratio(steel_area_val: float, concrete_area_val: float) -> float:
    """
    Return the ratio of steel area to concrete area.

    Args:
        steel_area_val (float): Total steel area (must be positive).
        concrete_area_val (float): Concrete area (must be positive).

    Returns:
        float: The ratio of steel to concrete area.

    Raises:
        ValueError: If 'steel_area_val' <= 0 or 'concrete_area_val' <= 0.
    """
    if steel_area_val <= 0:
        raise ValueError("Steel area must be positive.")
    if concrete_area_val <= 0:
        raise ValueError("Concrete area must be positive.")
    return steel_area_val / concrete_area_val


def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results


def transverseSpacing(length,db,ds,nbar,segment,cover=40):
    '''
    All units are in milimeters. 
    bc = 
    db = 
    ds = 
    nbar = 
    segment = no. of segment 
    '''
    wyi = (length - (cover*2) -( ds*2) - (db * nbar)) / segment  # fiber strip thickness along y
    return wyi

def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results

def clearTransverSpace(length,n,db,dbs,cover,segment):
    value_init = length - (2*cover) - (db*n) - (dbs*2) 
    value = value_init/segment
    return value

# #The Effective Lateral Confining Stress Concrete 
def compute_transverse_area(nlegs, db):
    """
    Compute the total cross-sectional area of a single set of ties/hoops,
    given the number of legs and the bar diameter (both in millimeters).

    Parameters
    ----------
    nlegs : int
        Number of bar legs in the tie/hoop set.
    bar_diam_mm : float
        Diameter of each bar (mm).

    Returns
    -------
    float
        Total cross-sectional area in mm^2.
    """
    # Area of one bar in mm^2
    area_bars = (math.pi / 4.0) * math.pow(db,2)
    # Multiply by number of legs
    total_area = nlegs * area_bars
    return total_area

def compute_ratio_trans(area_steel,core_length,spacing):
    results = area_steel/(core_length*spacing)
    return results

def compute_effective_conf_stress(ke,rho,fyh):
    fl = ke*rho*fyh 
    return fl

class generateFIBERMGT:
    '''
    name = name of element 
    fc = compressive strength of concrete (MPa)
    length_x = length of section along x 
    length_y = length of section along y
    db = main reinforcement diameter
    ds = secondary reinforcement diameter 
    ndbx = number of bars along x
    ndby = number of bars along y 
    nsegx = number of segment along x 
    nsegy = number of segment along y 
    cover = 40 , in millimeters  
    '''
    def __init__(self,name,fc,fy,length_x,length_y,total_n,db,ds,ndbx,ndby,nsegx,nsegy,spacing,spacing_prime,cover):
        self.name = name 
        self.fc= fc 
        self.fy= fy 
        self.length_x= length_x #core dimension at x
        self.length_y= length_y #core dimension at y   
        self.total_n = total_n          
        self.db = db 
        self.ds = ds
        self.ndbx= ndbx # No. of bars at X
        self.ndby= ndby # No. of bars a Y
        self.nsegx= nsegx #
        self.nsegy= nsegy #
        self.spacing= spacing #
        self.spacing_prime= spacing_prime #        
        self.cover = cover

    def unconfinedConcreteData(self):

        MPATOKPA = 1000
        concrete_expect = 1.5 
        fco_prime = concrete_expect*self.fc*MPATOKPA #kN/m      
        self.fco_prime = fco_prime
        # print(fco_prime)
        eco = 0.002 # Unconfined Concrete Strain (default)
        ecy = 0.0014 #Yield Strain for Unconfined Concrete (default) 
        esp = 0.02 # Spalling Stain for Unconfined Concrete (default)

        
        
        ec = 5000*math.sqrt(fco_prime/1000)*MPATOKPA #Elastic Modulus of Concrete 
        ft = 0.62*math.sqrt(fco_prime/1000)*MPATOKPA #Tensile Strength concrete 
        et = ft/ec
        self.eco = eco
        self.ecy = ecy
        self.esp = esp
        self.ec = ec
        self.ft = ft
        self.et = et
        data = f'{self.name},CONC, MANDER, 1, YES, {fco_prime}, NO, {self.eco}, NO, {self.ecy}, NO, {self.esp}, 0, {self.ec}, 0, {self.ft}, 0, {self.et},'
        print(data)
        return data
    def sectiondata(self):
        #Concrete Core Dimension to Center Line of Perimeter hoops
       
       
        core_x = concrete_core(self.length_x,self.cover,self.ds)
        core_y = concrete_core(self.length_y,self.cover,self.ds)
        self.core_x = core_x
        self.core_y = core_y

        print("bc = ", core_x)
        print("dc = ", core_y)
        #Concrete Core Dimension to Center Line of Perimeter hoops
        wxi = transverseSpacing(self.length_x,self.db,self.ds,self.ndbx,self.nsegx,self.cover)
        wyi = transverseSpacing(self.length_y,self.db,self.ds,self.ndby,self.nsegy,self.cover)
        print("wxi = ", wxi)
        print("wyi = ", wyi)
        seg_x_arr = []
        seg_y_arr = []
        
        n_seg_x_length = self.nsegx 
        n_seg_y_length = self.nsegy 

        if(n_seg_x_length == 0):
            n_seg_x_length = 1
        if(n_seg_y_length == 0):
            n_seg_y_length = 1
        
        length_x = n_seg_x_length
        length_y = n_seg_y_length

        i = 0 
        while i < length_x:
            # print(i)
            seg_x_arr.append(wxi)
            i += 1
        j = 0 
        while j < length_y: 
            # print(j)
            seg_y_arr.append(wyi)
            j += 1
        self.seg_x_arr_output = seg_x_arr
        self.seg_y_arr_output = seg_y_arr
        print(seg_x_arr)
        print(seg_y_arr)
    def confined_core_areas(self):
        """
        Compute:
        1) Ac  -- Core area
        2) Acc -- Effective concrete core area (Ac * (1 - rho_cc))
        3) Ae  -- Effectively confined core area
        4) kg  -- Ae / Ac
        5) ke  -- Ae / Acc
        
        Parameters
        ----------
        bc : float
            Core dimension in one direction (centerline of hoops). (mm)
        dc : float
            Core dimension in the orthogonal direction. (mm)
        rho_cc : float
            Longitudinal reinforcement ratio over the core (dimensionless).
        s_prime : float
            Clear spacing inside hoops or ties (same units as bc, dc).
        w_y : list of floats
            Clear tie widths in the y-direction (the w'_y_i values). (mm)
        w_z : list of floats
            Clear tie widths in the z-direction (the w'_z_j values). (mm)
        
        Returns
        -------
        Ac : float
            Core area = bc * dc
        Acc : float
            Effective concrete core area = Ac * (1 - rho_cc)
        Ae : float
            Effectively confined core area
        kg : float
            Confinement effectiveness coeff. = Ae / Ac
        ke : float
            Alternative confinement coeff.   = Ae / Acc
        """
        MPATOKPA = 1000
        fy_expected = self.fy*1.25*MPATOKPA
        self.fy_expected = fy_expected
        TOM2 = 1/(1000*1000)        
        #Transverse area
        trans_area_x = compute_transverse_area(self.ndbx, self.ds)
        trans_area_y = compute_transverse_area(self.ndby, self.ds)

        psx = compute_ratio_trans(trans_area_x,self.core_x,self.spacing)
        psy = compute_ratio_trans(trans_area_y,self.core_y,self.spacing)

        # 1) Core area
        area_core = self.core_x * self.core_y
        self.core = area_core*TOM2
        # print(Ac)
        dbm_area = area_diam(self.db)
        total_area = steel_area(self.total_n,dbm_area)
        rho_core = area_ratio(total_area,area_core)
        print("core ratio : ",rho_core)

        dbs_area = area_diam(self.ds)*TOM2 #MGT Report 
        self.dbs_area = dbs_area
        print("ds area : ",dbs_area)
        
        
        #Total Area of Confinement Rebars 
        area_conf_x = dbs_area*self.ndbx #MGT Report
        area_conf_y = dbs_area*self.ndby #MGT Report

        print("conf x-dir : ",area_conf_x)
        print("conf x-dir : ",area_conf_y)
        # 2) Effective concrete core area
        acc = area_core * (1.0 - rho_core)
        


        # Sums of squared tie widths
        sum_wy2 = sum(math.pow(wi,2)/6 for wi in self.seg_x_arr_output)
        sum_wz2 = sum(math.pow(wj,2)/6 for wj in self.seg_y_arr_output)  
        # 3) Ae (effectively confined core area)
        #    Matches the bracketed expression:
        #    2*( (sum wy_i^2)/6 + (sum wz_j^2)/6 ) ...
        factor = 2*(sum_wy2 + sum_wz2)
        Ae = (area_core - factor) * ((1.0 - self.spacing_prime/(2.0*self.core_x)) * (1.0 - self.spacing_prime/(2.0*self.core_y)))
        
        # 4) kg = Ae / Ac
        kg = Ae / area_core
        

        ke = Ae / acc if acc != 0.0 else float('nan')
        
        #The Effective Lateral Confining Stress on the Concrete 
        flx = compute_effective_conf_stress(ke,psx,fy_expected) #MGT REPORT
        fly = compute_effective_conf_stress(ke,psy,fy_expected) #MGT REPORT
        self.flx = flx
        self.fly = fly       
        print("flx : ",flx) 
        print("fly : ",fly)


        results = {
            "ac" : area_core,
            "acc" : acc*TOM2, #m2
            "ae" : Ae*TOM2, #m2
            "kg" : kg,
            "ke" : ke 
        }

        # self.ac = acc*TOM2
        self.acc = acc*TOM2
        self.ae = Ae*TOM2
        self.kg = kg
        self.ke = ke    
        print("ac",area_core)
        print("acc",acc*TOM2)
        print("ae",Ae*TOM2)
        print("kg",kg)
        print("ke",ke)

        return results
    def confined_concrete_strength_and_strain(self):
        """
        Compute the confined concrete strength (f'cc) and strain (eps_cc)
        using the six-step procedure shown in your reference.

        Steps:
        1) q = f'l1 / f'l2  (with f'l2 >= f'l1)
        2) A = 6.886 - [ (0.6069 + 17.275*q ) * exp(-4.989*q) ]
        3) B = (4.5 / 5)*[ 0.9849 - 0.6306*exp(-3.8939*q) ]^(-5) - 0.1
            -- (exact exponent form may vary; confirm with your text!)
        4) x' = ( f'l1 + f'l2 ) / ( 2 * f'co )
        5) k1 = A * [ 0.1 + 0.9 / (1 + B * x') ]
        6) f'cc = f'co * [ 1 + k1 * x' ]
            eps_cc = eps_co * [ 1 + 5 * ( f'cc / f'co - 1 ) ]

        Parameters
        ----------
        fl1 : float
            f'l1, smaller lateral confining stress (kN/m²).
        fl2 : float
            f'l2, larger lateral confining stress (kN/m²). Must satisfy fl2 >= fl1.
        fco : float
            Unconfined concrete strength (kN/m²).
        eco : float
            Unconfined concrete strain at f'co (dimensionless, ~0.002 typical).

        Returns
        -------
        fcc : float
            Confined concrete strength, f'cc (kN/m²).
        ecc : float
            Strain at f'cc (dimensionless).
        """
        # MPATOKPA = 1000
        # concrete_expect = 1.5 
        # fco_prime = concrete_expect*self.fc*MPATOKPA #kN/m     

        eco = 0.002
        # 1) q = f'l1 / f'l2 (assuming fl2 >= fl1)
        #    If fl2 < fl1 in your data, swap them or check:
        q = self.flx / self.fly if self.fly != 0 else 0.0

        # 2) A
        A = 6.886 - (0.6069 + 17.275*q) * math.exp(-4.989*q)

        # 3) B  (this expression is inferred from your snippet;
        #        please confirm exact exponent, etc., from your reference.)
        B = (4.5 / 5.0)*(0.9849 - 0.6306*math.exp(-3.8939*q))**(-5) - 0.1

        # 4) x' = (f'l1 + f'l2) / (2 * fco)
        x_prime = (self.flx + self.fly) / (2.0 * self.fco_prime)

        # 5) k1 = A * [ 0.1 + 0.9 / (1 + B * x') ]
        k1 = A * (0.1 + 0.9/(1.0 + B*x_prime))

        # 6) f'cc = f'co [1 + k1 x']
        fcc = self.fco_prime * (1.0 + k1*x_prime)

        #    eps_cc = eps_co [1 + 5 (f'cc / f'co - 1)]
        ecc = eco * (1.0 + 5.0*(fcc/self.fco_prime - 1.0))
        print("fcc = ",fcc)
        print("ecc = ",ecc)
        self.fcc = fcc 
        self.ecc = ecc
        # return fcc, ecc
    # def generateMGT(self):
        #Concrete Type
        #Confinement Effectivenes Coefficeint 
        #
        # result1 = f'   ${self.name}, CONC, MANDER, 1, YES, ${self.fco_prime}, NO, {self.eco}, NO, {self.ecy}, NO, {self.esp}, 0, {self.ec}, {self.ft}, {self.et},'
        # result2 = f'0, NO, {self.core}, NO, {self.ae}, NO, {self.ke}, NO, 0, NO, {self.fy_expected},'
        # result3 = f'NO, 0, NO, {self.fcc}, NO, {self.ecc}, NO, {self.ecc}, NO, NO, 0, 0, NO, 0, 0, {}, {}, {}, {}, {}, {}'        
        #   C3000_C1_CONF, CONC, MANDER, 1, YES, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0795894, NO, 0.0205812, NO, 0.258592, NO, 0, NO, 28359, NO, 0.00296794, NO, 0.00207756, NO, NO, 0, 0, NO, 0, 1, 1, 0.16, 0.118, 1, 0.51, 0.105, 4, 0, 0, 10, D16, 0.0020106, NO, 0.0246397, 0, D10, 7.854e-05, 0.2, 0.19, 5, 2, 0.0003927, 0.00015708, NO, 344750, NO, NO, NO, 0.00385, NO, 0.00490875, NO, 343.227, NO, 437.614, NO, 0, NO, NO, 0, PUSHOVER



test = generateFIBERMGT("CHECK",20.68,276,250,600,10,16,10,2,5,1,4,200,190,40)
test.unconfinedConcreteData()
test.sectiondata()
test.confined_core_areas()
test.confined_concrete_strength_and_strain()

# '''
#   NO, 0, NO, 33534.9, NO, 0.00281074, NO, 0.00281074, NO, NO, 0, 0, NO, 0, 1, 1, 0.16, 0.118, 1, 0.51, 0.105, 4, 0, 0, 10, D16, 0.0020106, NO, 0.0246397, 0, D10, 7.854e-05, 0.2, 0.19, 5, 2, 0.0003927, 0.00015708, NO, 344750, NO, NO, NO, 0.00385, NO, 0.00490875, NO, 343.227, NO, 437.614, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C2_CONF, CONC, MANDER, 1, YES, 31020, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.78478e+07, 0, 3453.13, 0.000124, 0, NO, 0.0434915, NO, 0.00324932, NO, 0.0747116, NO, 0, NO, 31889.4, NO, 0.00228026, NO, 0.00228026, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.112, 3, 0, 0, 8, D16, 0.00160848, NO, 0.0356647, 0, D10, 7.854e-05, 0.2, 0.19, 4, 2, 0.00031416, 0.00015708, NO, 344750, NO, NO, NO, 0.00383122, NO, 0.00714, NO, 98.68, NO, 183.904, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C3_CONF, CONC, MANDER, 1, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0426873, NO, 0.00375611, NO, 0.0879914, NO, 0, NO, 27131.8, NO, 0.00249331, NO, 0.00174531, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.068, 5, 0, 0, 12, D16, 0.00241272, NO, 0.0534971, 0, D10, 7.854e-05, 0.2, 0.19, 6, 2, 0.00047124, 0.00015708, NO, 344750, NO, NO, NO, 0.00574683, NO, 0.00714, NO, 174.33, NO, 216.592, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C1_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C2_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C3_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   GRADE 40 EXPECTED, STEEL, PM, 344750, 475000, 200000000, 0.015, 0.12, PUSHOVER


# '''