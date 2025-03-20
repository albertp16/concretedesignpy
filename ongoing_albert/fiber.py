# unconcrete 
import math
# import concretedesignpy
fc = 20.68 #Mpa
MPATOKPA = 1000
concrete_expect = 1.5 





fco_prime = concrete_expect*fc*MPATOKPA #kN/m

print(fco_prime)

'''
  C3000_C1_CONF, CONC, MANDER, 1, YES, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0795894, NO, 0.0205812, NO, 0.258592, NO, 0, NO, 28359, NO, 0.00296794, NO, 0.00207756, NO, NO, 0, 0, NO, 0, 1, 1, 0.16, 0.118, 1, 0.51, 0.105, 4, 0, 0, 10, D16, 0.0020106, NO, 0.0246397, 0, D10, 7.854e-05, 0.2, 0.19, 5, 2, 0.0003927, 0.00015708, NO, 344750, NO, NO, NO, 0.00385, NO, 0.00490875, NO, 343.227, NO, 437.614, NO, 0, NO, NO, 0, PUSHOVER
  C3000_C2_CONF, CONC, MANDER, 1, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0434915, NO, 0.00324932, NO, 0.0747116, NO, 0, NO, 26724.1, NO, 0.00233565, NO, 0.00163495, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.112, 3, 0, 0, 8, D16, 0.00160848, NO, 0.0356647, 0, D10, 7.854e-05, 0.2, 0.19, 4, 2, 0.00031416, 0.00015708, NO, 344750, NO, NO, NO, 0.00383122, NO, 0.00714, NO, 98.68, NO, 183.904, NO, 0, NO, NO, 0, PUSHOVER
  C3000_C3_CONF, CONC, MANDER, 1, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0.0426873, NO, 0.00375611, NO, 0.0879914, NO, 0, NO, 27131.8, NO, 0.00249331, NO, 0.00174531, NO, NO, 0, 0, NO, 0, 1, 1, 0.11, 0.068, 1, 0.41, 0.068, 5, 0, 0, 12, D16, 0.00241272, NO, 0.0534971, 0, D10, 7.854e-05, 0.2, 0.19, 6, 2, 0.00047124, 0.00015708, NO, 344750, NO, NO, NO, 0.00574683, NO, 0.00714, NO, 174.33, NO, 216.592, NO, 0, NO, NO, 0, PUSHOVER
  C3000_C1_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
  C3000_C2_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
  C3000_C3_UNCONF, CONC, MANDER, 0, NO, 25856.3, NO, 0.002, NO, 0.0014, NO, 0.02, 0, 2.54245e+07, 0, 3152.64, 0.000124, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER
  GRADE 40 EXPECTED, STEEL, PM, 344750, 475000, 200000000, 0.015, 0.12, PUSHOVER

'''
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


db = 16 
ds = 10
n = 10
cx = 250
cy = 600 
cover = 40
section_area = cx*cy

def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results


#obtaining the Concrete Core Dimenions 
bcx = concrete_core(cx,cover,ds)
bcy = concrete_core(cy,cover,ds)
# print(bcx)
# print(bcy)
section_area_core = bcx*bcy
print("section_area_core : " + str(section_area_core))
#Transver spacing 
def transverseSpacing(bc,db,ds,nbar,segment):
    wyi = (bc - ds - (db * nbar)) / segment  # fiber strip thickness along y
    return wyi

wxi = transverseSpacing(bcx,db,ds,2,1)
wyi = transverseSpacing(bcy,db,ds,5,4)
print("wxi : " + str(wxi))
print("wyi : " + str(wyi))


db16 = area_diam(db)
total_area = steel_area(n,db16)
rho_core = area_ratio(total_area,section_area_core)
print(total_area)
print(rho_core)

acc = section_area_core * (1.0 - rho_core) #TODO to convert ot m2
print("acc : " + str(acc)) 

def confined_core_areas(bc, dc, rho_cc, s_prime, w_y, w_z):
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
        Core dimension in one direction (centerline of hoops).
    dc : float
        Core dimension in the orthogonal direction.
    rho_cc : float
        Longitudinal reinforcement ratio over the core (dimensionless).
    s_prime : float
        Clear spacing inside hoops or ties (same units as bc, dc).
    w_y : list of floats
        Clear tie widths in the y-direction (the w'_y_i values).
    w_z : list of floats
        Clear tie widths in the z-direction (the w'_z_j values).
    
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
    
    # 1) Core area
    Ac = bc * dc
    # print(Ac)
    # 2) Effective concrete core area
    Acc = Ac * (1.0 - rho_cc)
    
    # Sums of squared tie widths
    sum_wy2 = sum(math.pow(wi,2)/6 for wi in w_y)
    sum_wz2 = sum(math.pow(wj,2)/6 for wj in w_z)
    
    # 3) Ae (effectively confined core area)
    #    Matches the bracketed expression:
    #    2*( (sum wy_i^2)/6 + (sum wz_j^2)/6 ) ...
    factor = 2*(sum_wy2 + sum_wz2)
    Ae = (Ac - factor) * ((1.0 - s_prime/(2.0*bc)) * (1.0 - s_prime/(2.0*dc)))
    
    # 4) kg = Ae / Ac
    kg = Ae / Ac
    
    #    ke = Ae / Acc = Ae / (Ac * (1 - rho_cc))
    ke = Ae / Acc if Acc != 0.0 else float('nan')
    results = {
        "ac" : Ac,
        "acc" : Acc,
        "ae" : Ae,
        "kg" : kg,
        "ke" : ke 
    }
    return results


core_data = confined_core_areas(bcx, bcy, rho_core, 190, [wxi], [wyi,wyi,wyi,wyi])
print("ac  -> " ,core_data["ac"])
print("acc -> " ,core_data["acc"])
print("ae  -> " ,core_data["ae"])
print("kg  -> " ,core_data["kg"])
print("ke  -> " ,core_data["ke"])


def concrete_core(length,cover,dia_stirups):
    results = length - cover - cover - (dia_stirups/2) - (dia_stirups/2)
    return results

def clearTransverSpace(length,n,db,dbs,cover,segment):
    value_init = length - (2*cover) - (db*n) - (dbs*2) 
    value = value_init/segment
    return value

def generateMGTCONC(name,fc,fy,type,cx,cy,ds,cover,cxn,cyn):
    # fc = 20.68 #Mpa
    MPATOKPA = 1000
    concrete_expect = 1.5 
    fco_prime = concrete_expect*fc*MPATOKPA #kN/
    
    core_x = cx - cover - cover - (ds/2) - (ds/2)
    core_y = cx - cover - cover - (ds/2) - (ds/2)


     # print(fco_prime)
    eco = 0.002 # Unconfined Concrete Strain
    ecy = 0.0014 #Yield Strain for Unconfined Concrete 
    esp = 0.02 # Spalling Stain for Unconfined Concrete

    ec = 5000*math.sqrt(fco_prime/1000)*MPATOKPA #Elastic Modulus of Concrete 
    ft = 0.62*math.sqrt(fco_prime/1000)*MPATOKPA #Tensile Strength concrete
    et = ft/ec
    if type == True:
        data = f'{name},CONC, MANDER, 1, YES, {fco_prime}, NO, {eco}, NO, {ecy}, NO, {esp}, 0, {ec}, 0, {ft}, 0, {et},'
    else:
        data = f'{name},CONC, MANDER, 1, YES, {fco_prime}, NO, {eco}, NO, {ecy}, NO, {esp}, 0, {ec}, 0, {ft}, 0, {et},  0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, 0, NO, 0, 1, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, , 0, NO, 0, 0, , 0, 0, 0, 0, 0, 0, 0, NO, 0, NO, NO, NO, 0, NO, 0, NO, 0, NO, 0, NO, 0, NO, NO, 0, PUSHOVER'        
    print(data)

# generateMGTCONC("C3000_C1_CONF",fc,276,True)
# generateMGTCONC("C3000_UNCONF",fc,276,True)