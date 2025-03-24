import math
import pandas as pd

def eRF(exposure,type):
    data = {
        "interior" : {
            "carbon" : 0.95,
            "glass" : 0.75,
            "aramid" : 0.85
        },
        "exterior" : {
            "carbon" : 0.85,
            "glass" : 0.65,
            "aramid" : 0.75
        },
        "aggressive" : {
            "carbon" : 0.85,
            "glass" : 0.50,
            "aramid" : 0.70
            },
    }
    results = data[exposure][type]
    return results

def beta_1(fc): 
    ##NSCP 2015 422.2.2.4.3
    value = 1.05 - 0.05*(fc/1000)
    if(17 <= fc and fc <= 28):
        value = 0.85
    elif (28 <= fc and fc < 55):
        value = 0.85 - ((0.05*(fc-28))/7)
    else:
        value = 0.65

    return value
    
def areaCircle(db):
    pi = math.pi
    value = (pi/4)*(math.pow(db,2))
    return value

ebi = ( MDL * (df - (k*d)) )  /( Icr * ec_convert ) 

def cSolver():## //Numerical Solver
        
    ##--> Data for the tabulations
    epsilon_fe_arr = []
    espilon_c_arr = []
    eespilon_s_arr = []
    e_c_prime_arr = []
    beta_one_arr = []
    alpha_one_arr = []
    c_init_arr = []

    c_init = 0.20*input['d'] ##initialize depth
    i = 0
    
    while i <= 20:
       
        epsilon_fe = 0.003*( (input['h'] - c_init) /c_init) - ebi
        espilon_c = (0.009 + ebi)*(c_init/(input['h'] - c_init))
        espilon_s = (0.009 + ebi)*((input['d'] - c_init)/(input['h'] - c_init))

        fs = input['Es']*espilon_s
        fs = 414 ##TODO fix on the yields
            
        f_fe = input['frp']['Ef']*0.009
            
        e_c_prime = 1.7*input['fc']/ec
        beta_one = ((4*e_c_prime )- espilon_c)/((6*e_c_prime)-(2*espilon_c))
        alpha_one = ((3*e_c_prime*espilon_c) - math.pow(espilon_c,2))/(3*beta_one*math.pow(e_c_prime,2))    
        c_init = (area_steel*fs + af*f_fe)/(alpha_one*input['fc']*beta_one*input['w'])

        epsilon_fe_arr.append(epsilon_fe)
        espilon_c_arr.append(espilon_c)
        eespilon_s_arr.append(espilon_s)
        e_c_prime_arr.append(e_c_prime)
        beta_one_arr.append(beta_one)
        alpha_one_arr.append(alpha_one)
        c_init_arr.append(c_init)
    
        c = c_init
        i += 1
    results = {
        "value" : c,
        "epsilon_fe" : epsilon_fe_arr,
        "espilon_c" : espilon_c_arr,
        "espilon_s" : eespilon_s_arr,
        "e_c_prime" : e_c_prime_arr,
        "beta_one" : beta_one_arr,
        "alpha_one" : alpha_one_arr,
        "c_init" : c_init_arr,
    }
        
    return results

c_data = cSolver()
df = pd.DataFrame({
      'c (mm)': c_data['c_init'],
      'ε_fe': c_data['epsilon_fe'],    
      'ε_c': c_data['espilon_c'],
      'ε_s': c_data['espilon_s'],
      'ε_c\'': c_data['e_c_prime'],
      'β1': c_data['beta_one'],
      'α': c_data['alpha_one'],
      })