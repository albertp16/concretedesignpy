# unconcrete 
import math
# # import concretedesignpy
# fc = 20.68 #Mpa
# fy = 276
# MPATOKPA = 1000
# concrete_expect = 1.5 
# steel_expect = 1.25 




# fco_prime = concrete_expect*fc*MPATOKPA #kN/m
# fyh = steel_expect*fy*MPATOKPA #kN/m
# print("fco -> ",fco_prime)
# print("fyh -> ",fyh)

# '''
#   C3000_C1_CONF, CONC, MANDER, 1, YES, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0795894, NO, 0.0205812, NO, 0.258592, NO, 0, NO, 28359, NO, 0.00296794, NO, 0.00207756, NO, NO, 0, 0, NO, 0, 1, 1, 0.16, 0.118, 1, 0.51, 0.105, 4, 0, 0, 10, D16, 0.0020106, NO, 0.0246397, 0, D10, 7.854e-05, 0.2, 0.19, 5, 2, 0.0003927, 0.00015708, NO, 344750, NO, NO, NO, 0.00385, NO, 0.00490875, NO, 343.227, NO, 437.614, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C2_CONF, CONC, MANDER, 1, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0434915, NO, 0.00324932, NO, 0.0747116, NO, 0, NO, 26724.1, NO, 0.00233565, NO, 0.00163495, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.112, 3, 0, 0, 8, D16, 0.00160848, NO, 0.0356647, 0, D10, 7.854e-05, 0.2, 0.19, 4, 2, 0.00031416, 0.00015708, NO, 344750, NO, NO, NO, 0.00383122, NO, 0.00714, NO, 98.68, NO, 183.904, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C3_CONF, CONC, MANDER, 1, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0426873, NO, 0.00375611, NO, 0.0879914, NO, 0, NO, 27131.8, NO, 0.00249331, NO, 0.00174531, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.068, 5, 0, 0, 12, D16, 0.00241272, NO, 0.0534971, 0, D10, 7.854e-05, 0.2, 0.19, 6, 2, 0.00047124, 0.00015708, NO, 344750, NO, NO, NO, 0.00574683, NO, 0.00714, NO, 174.33, NO, 216.592, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C1_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C2_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   C3000_C3_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
#   GRADE 40 EXPECTED, STEEL, PM, 344750, 475000, 200000000, 0.015, 0.12, PUSHOVER

# '''
# import math

# def area_diam(diameter: float) -> float:
#     """
#     Return the cross-sectional area of a circle given its diameter.

#     Args:
#         diameter (float): The diameter of the circle (must be positive).

#     Returns:
#         float: The area of the circle.

#     Raises:
#         ValueError: If 'diameter' is not positive.
#     """
#     if diameter <= 0:
#         raise ValueError("Diameter must be positive.")
#     return (math.pi / 4) * (diameter ** 2)


# def steel_area(num_bars: int, bar_area: float) -> float:
#     """
#     Return total steel area for a given number of bars.

#     Args:
#         num_bars (int): Number of bars (must be positive).
#         bar_area (float): Cross-sectional area of one bar (must be positive).

#     Returns:
#         float: The total steel area.

#     Raises:
#         ValueError: If 'num_bars' <= 0 or 'bar_area' <= 0.
#     """
#     if num_bars <= 0:
#         raise ValueError("Number of bars must be positive.")
#     if bar_area <= 0:
#         raise ValueError("Bar area must be positive.")
#     return num_bars * bar_area


# def area_ratio(steel_area_val: float, concrete_area_val: float) -> float:
#     """
#     Return the ratio of steel area to concrete area.

#     Args:
#         steel_area_val (float): Total steel area (must be positive).
#         concrete_area_val (float): Concrete area (must be positive).

#     Returns:
#         float: The ratio of steel to concrete area.

#     Raises:
#         ValueError: If 'steel_area_val' <= 0 or 'concrete_area_val' <= 0.
#     """
#     if steel_area_val <= 0:
#         raise ValueError("Steel area must be positive.")
#     if concrete_area_val <= 0:
#         raise ValueError("Concrete area must be positive.")
#     return steel_area_val / concrete_area_val


# db = 16 
# ds = 10
# n = 10
# cx = 250
# cy = 600 
# cover = 40
# section_area = cx*cy

def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results


# #obtaining the Concrete Core Dimenions 
# bcx = concrete_core(cx,cover,ds)
# bcy = concrete_core(cy,cover,ds)
# # print(bcx)
# # print(bcy)
# section_area_core = bcx*bcy
# print("section_area_core : " + str(section_area_core))
# #Transver spacing 
def transverseSpacing(bc,db,ds,nbar,segment):
    wyi = (bc - ds - (db * nbar)) / segment  # fiber strip thickness along y
    return wyi

# wxi = transverseSpacing(bcx,db,ds,2,1)
# wyi = transverseSpacing(bcy,db,ds,5,4)
# print("wxi : " + str(wxi))
# print("wyi : " + str(wyi))


# db16 = area_diam(db)
# total_area = steel_area(n,db16)
# rho_core = area_ratio(total_area,section_area_core)
# print(total_area)
# print(rho_core)

# acc = section_area_core * (1.0 - rho_core) #TODO to convert ot m2
# print("acc : " + str(acc)) 

# def confined_core_areas(bc, dc, rho_cc, s_prime, w_y, w_z):
#     """
#     Compute:
#       1) Ac  -- Core area
#       2) Acc -- Effective concrete core area (Ac * (1 - rho_cc))
#       3) Ae  -- Effectively confined core area
#       4) kg  -- Ae / Ac
#       5) ke  -- Ae / Acc
    
#     Parameters
#     ----------
#     bc : float
#         Core dimension in one direction (centerline of hoops). (mm)
#     dc : float
#         Core dimension in the orthogonal direction. (mm)
#     rho_cc : float
#         Longitudinal reinforcement ratio over the core (dimensionless).
#     s_prime : float
#         Clear spacing inside hoops or ties (same units as bc, dc).
#     w_y : list of floats
#         Clear tie widths in the y-direction (the w'_y_i values). (mm)
#     w_z : list of floats
#         Clear tie widths in the z-direction (the w'_z_j values). (mm)
    
#     Returns
#     -------
#     Ac : float
#         Core area = bc * dc
#     Acc : float
#         Effective concrete core area = Ac * (1 - rho_cc)
#     Ae : float
#         Effectively confined core area
#     kg : float
#         Confinement effectiveness coeff. = Ae / Ac
#     ke : float
#         Alternative confinement coeff.   = Ae / Acc
#     """
    
#     # 1) Core area
#     Ac = bc * dc
#     # print(Ac)
#     # 2) Effective concrete core area
#     Acc = Ac * (1.0 - rho_cc)
    
#     # Sums of squared tie widths
#     sum_wy2 = sum(math.pow(wi,2)/6 for wi in w_y)
#     sum_wz2 = sum(math.pow(wj,2)/6 for wj in w_z)
    
#     # 3) Ae (effectively confined core area)
#     #    Matches the bracketed expression:
#     #    2*( (sum wy_i^2)/6 + (sum wz_j^2)/6 ) ...
#     factor = 2*(sum_wy2 + sum_wz2)
#     Ae = (Ac - factor) * ((1.0 - s_prime/(2.0*bc)) * (1.0 - s_prime/(2.0*dc)))
    
#     # 4) kg = Ae / Ac
#     kg = Ae / Ac
    
#     #    ke = Ae / Acc = Ae / (Ac * (1 - rho_cc))
#     ke = Ae / Acc if Acc != 0.0 else float('nan')
#     results = {
#         "ac" : Ac,
#         "acc" : Acc,
#         "ae" : Ae,
#         "kg" : kg,
#         "ke" : ke 
#     }
#     return results


# # core_data = confined_core_areas(bcx, bcy, rho_core, 190, [wxi], [wyi,wyi,wyi,wyi])
# # print("ac  -> " ,core_data["ac"])
# # print("acc -> " ,core_data["acc"])
# # print("ae  -> " ,core_data["ae"])
# # print("kg  -> " ,core_data["kg"])
# # print("ke  -> " ,core_data["ke"])


def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results

def clearTransverSpace(length,n,db,dbs,cover,segment):
    value_init = length - (2*cover) - (db*n) - (dbs*2) 
    value = value_init/segment
    return value

# #The Effective Lateral Confining Stress Concrete 
# def compute_transverse_area(nlegs, db):
#     """
#     Compute the total cross-sectional area of a single set of ties/hoops,
#     given the number of legs and the bar diameter (both in millimeters).

#     Parameters
#     ----------
#     nlegs : int
#         Number of bar legs in the tie/hoop set.
#     bar_diam_mm : float
#         Diameter of each bar (mm).

#     Returns
#     -------
#     float
#         Total cross-sectional area in mm^2.
#     """
#     # Area of one bar in mm^2
#     area_bars = (math.pi / 4.0) * math.pow(db,2)
#     # Multiply by number of legs
#     total_area = nlegs * area_bars
#     return total_area

# fyh_as_x = compute_transverse_area(2, ds)
# fyh_as_y = compute_transverse_area(5, ds)
# print(fyh_as_x)
# print(fyh_as_y)

# def compute_ratio_trans(area_steel,core_length,spacing):
#     results = area_steel/(core_length*spacing)
#     return results

# psx = compute_ratio_trans(fyh_as_x,bcx,200)
# psy = compute_ratio_trans(fyh_as_y,bcy,200)
# print(psx)
# print(psy)

# def compute_effective_conf_stress(ke,rho,fyh):
#     fl = ke*rho*fyh 
#     return fl

# def confined_concrete_strength_and_strain(
#     fl1, fl2,   # f'l1, f'l2 (lateral confining stresses in two directions)
#     fco,        # Unconfined concrete strength, f'co
#     eco=0.002         # Unconfined concrete strain at f'co (often ~0.002)
# ):
#     """
#     Compute the confined concrete strength (f'cc) and strain (eps_cc)
#     using the six-step procedure shown in your reference.

#     Steps:
#       1) q = f'l1 / f'l2  (with f'l2 >= f'l1)
#       2) A = 6.886 - [ (0.6069 + 17.275*q ) * exp(-4.989*q) ]
#       3) B = (4.5 / 5)*[ 0.9849 - 0.6306*exp(-3.8939*q) ]^(-5) - 0.1
#          -- (exact exponent form may vary; confirm with your text!)
#       4) x' = ( f'l1 + f'l2 ) / ( 2 * f'co )
#       5) k1 = A * [ 0.1 + 0.9 / (1 + B * x') ]
#       6) f'cc = f'co * [ 1 + k1 * x' ]
#          eps_cc = eps_co * [ 1 + 5 * ( f'cc / f'co - 1 ) ]

#     Parameters
#     ----------
#     fl1 : float
#         f'l1, smaller lateral confining stress (kN/m²).
#     fl2 : float
#         f'l2, larger lateral confining stress (kN/m²). Must satisfy fl2 >= fl1.
#     fco : float
#         Unconfined concrete strength (kN/m²).
#     eco : float
#         Unconfined concrete strain at f'co (dimensionless, ~0.002 typical).

#     Returns
#     -------
#     fcc : float
#         Confined concrete strength, f'cc (kN/m²).
#     ecc : float
#         Strain at f'cc (dimensionless).
#     """

#     # 1) q = f'l1 / f'l2 (assuming fl2 >= fl1)
#     #    If fl2 < fl1 in your data, swap them or check:
#     q = fl1 / fl2 if fl2 != 0 else 0.0

#     # 2) A
#     A = 6.886 - (0.6069 + 17.275*q) * math.exp(-4.989*q)

#     # 3) B  (this expression is inferred from your snippet;
#     #        please confirm exact exponent, etc., from your reference.)
#     B = (4.5 / 5.0)*(0.9849 - 0.6306*math.exp(-3.8939*q))**(-5) - 0.1

#     # 4) x' = (f'l1 + f'l2) / (2 * fco)
#     x_prime = (fl1 + fl2) / (2.0 * fco)

#     # 5) k1 = A * [ 0.1 + 0.9 / (1 + B * x') ]
#     k1 = A * (0.1 + 0.9/(1.0 + B*x_prime))

#     # 6) f'cc = f'co [1 + k1 x']
#     fcc = fco * (1.0 + k1*x_prime)

#     #    eps_cc = eps_co [1 + 5 (f'cc / f'co - 1)]
#     ecc = eco * (1.0 + 5.0*(fcc/fco - 1.0))
#     print("fcc = ",fcc)
#     print("ecc = ",ecc)
#     return fcc, ecc

# test = confined_concrete_strength_and_strain(343.2266,437.6139,fco_prime,eco=0.002)

# def generateMGTCONC(name,fc,fy,type,cx,cy,ds,cover,cxn,cyn):
#     # fc = 20.68 #Mpa
#     MPATOKPA = 1000
#     concrete_expect = 1.5 
#     fco_prime = concrete_expect*fc*MPATOKPA #kN/
    
#     core_x = cx - cover - cover - (ds/2) - (ds/2)
#     core_y = cx - cover - cover - (ds/2) - (ds/2)


#      # print(fco_prime)
#     eco = 0.002 # Unconfined Concrete Strain
#     ecy = 0.0014 #Yield Strain for Unconfined Concrete 
#     esp = 0.02 # Spalling Stain for Unconfined Concrete

#     ec = 5000*math.sqrt(fco_prime/1000)*MPATOKPA #Elastic Modulus of Concrete 
#     ft = 0.62*math.sqrt(fco_prime/1000)*MPATOKPA #Tensile Strength concrete
#     et = ft/ec
#     if type == True:
#         data = f'{name},CONC, MANDER, 1, YES, {fco_prime}, NO, {eco}, NO, {ecy}, NO, {esp}, 0, {ec}, 0, {ft}, 0, {et},'
#     else:
#         data = f'{name},CONC, MANDER, 1, YES, {fco_prime}, NO, {eco}, NO, {ecy}, NO, {esp}, 0, {ec}, 0, {ft}, 0, {et},  0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER'        
#     print(data)

# generateMGTCONC("C3000_C1_CONF",fc,276,True)
# generateMGTCONC("C3000_UNCONF",fc,276,True)

class generateFIBERMGT:
    def __init__(self,name,fc,length_x,length_y,db,ds,ndbx,ndby,nsegx,nsegy,cover):
        self.name = name 
        self.fc= fc 
        self.length_x= length_x #core dimension at x
        self.length_y= length_y #core dimension at y    
        self.db = db 
        self.ds = ds
        self.ndbx= ndbx #
        self.ndby= ndby #
        self.nsegx= nsegx #
        self.nsegy= nsegy #
        self.cover = cover
    def unconfinedConcreteData(self):
        MPATOKPA = 1000
        concrete_expect = 1.5 
        fco_prime = concrete_expect*self.fc*MPATOKPA #kN/m        
        # print(fco_prime)
        eco = 0.002 # Unconfined Concrete Strain (default)
        ecy = 0.0014 #Yield Strain for Unconfined Concrete (default) 
        esp = 0.02 # Spalling Stain for Unconfined Concrete (default)
        
        ec = 5000*math.sqrt(fco_prime/1000)*MPATOKPA #Elastic Modulus of Concrete 
        ft = 0.62*math.sqrt(fco_prime/1000)*MPATOKPA #Tensile Strength concrete 
        et = ft/ec
        data = f'{self.name},CONC, MANDER, 1, YES, {fco_prime}, NO, {eco}, NO, {ecy}, NO, {esp}, 0, {ec}, 0, {ft}, 0, {et},'
        print(data)
        return data
    def section_data(self):
        #Concrete Core Dimension to Center Line of Perimeter hoops
        core_x = concrete_core(self.length_x,self.cover,self.ds)
        core_y = concrete_core(self.length_y,self.cover,self.ds)
        #Concrete Core Dimension to Center Line of Perimeter hoops
        wxi = transverseSpacing(self.length_x,self.db,self.ds,self.ndbx,self.nsegx)
        wyi = transverseSpacing(self.length_y,self.db,self.ds,self.ndby,self.nsegy)
        seg_x_arr = []
        seg_y_arr = []
        
        n_seg_x_length = self.nsegx - 1
        n_seg_y_length = self.nsegy - 1

        if(n_seg_x_length == 0):
            n_seg_x_length = 1
        if(n_seg_y_length == 0):
            n_seg_y_length = 1
        
        length_x = n_seg_x_length
        length_y = n_seg_y_length

        i = 0 
        while i < length_x:
            print(i)
            seg_x_arr.append(wxi)
            i += 1
        j = 0 
        while j < length_y:
            print(j)
            seg_y_arr.append(wyi)
            j += 1
        print(seg_x_arr)
        print(seg_y_arr)

    # def confinementEffectiveCoefficient(self):
    #     """
    #     Compute:
    #     1) Ac  -- Core area
    #     2) Acc -- Effective concrete core area (Ac * (1 - rho_cc))
    #     3) Ae  -- Effectively confined core area
    #     4) kg  -- Ae / Ac
    #     5) ke  -- Ae / Acc
        
    #     Parameters
    #     ----------
    #     bc : float
    #         Core dimension in one direction (centerline of hoops). (mm)
    #     dc : float
    #         Core dimension in the orthogonal direction. (mm)
    #     rho_cc : float
    #         Longitudinal reinforcement ratio over the core (dimensionless).
    #     s_prime : float
    #         Clear spacing inside hoops or ties (same units as bc, dc).
    #     w_y : list of floats
    #         Clear tie widths in the y-direction (the w'_y_i values). (mm)
    #     w_z : list of floats
    #         Clear tie widths in the z-direction (the w'_z_j values). (mm)
        
    #     Returns
    #     -------
    #     Ac : float
    #         Core area = bc * dc
    #     Acc : float
    #         Effective concrete core area = Ac * (1 - rho_cc)
    #     Ae : float
    #         Effectively confined core area
    #     kg : float
    #         Confinement effectiveness coeff. = Ae / Ac
    #     ke : float
    #         Alternative confinement coeff.   = Ae / Acc
    #     """
        
    #     # 1) Core area
    #     area_core = self.bc * self.dc
    #     # print(Ac)
    #     # 2) Effective concrete core area
    #     Acc = area_core * (1.0 - rho_cc)
        
    #     # Sums of squared tie widths
    #     sum_wy2 = sum(math.pow(wi,2)/6 for wi in w_y)
    #     sum_wz2 = sum(math.pow(wj,2)/6 for wj in w_z)
        
    #     # 3) Ae (effectively confined core area)
    #     #    Matches the bracketed expression:
    #     #    2*( (sum wy_i^2)/6 + (sum wz_j^2)/6 ) ...
    #     factor = 2*(sum_wy2 + sum_wz2)
    #     Ae = (area_core - factor) * ((1.0 - s_prime/(2.0*bc)) * (1.0 - s_prime/(2.0*dc)))
        
    #     # 4) kg = Ae / Ac
    #     kg = Ae / area_core
        
    #     #    ke = Ae / Acc = Ae / (Ac * (1 - rho_cc))
    #     ke = Ae / Acc if Acc != 0.0 else float('nan')
    #     results = {
    #         "ac" : area_core,
    #         "acc" : Acc,
    #         "ae" : Ae,
    #         "kg" : kg,
    #         "ke" : ke 
    #     }
    #     return results        
    #     # return 0

test = generateFIBERMGT("CHECK",20.68,600,250,16,10,5,2,1,4,40)
test.unconfinedConcreteData()
test.section_data()