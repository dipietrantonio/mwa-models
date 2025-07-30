from models import dispersive_delay


def compute_parameters(centre_fq_mhz, dm, tres):
    
    cf_ghz = centre_fq_mhz / 1000
    bwhz = 0.03072
    f_lo_ghz = cf_ghz - bwhz / 2
    
    delta_ghz = 0.00001 # 10 kHz
    f_hi_ghz = f_lo_ghz + delta_ghz
    
    delay = dispersive_delay(dm, f_lo_ghz, f_hi_ghz) / 1000
    while delay < tres:
        f_hi_ghz += delta_ghz
        delay = dispersive_delay(dm, f_lo_ghz, f_hi_ghz) / 1000
    
    print("Final bandwidth required (kHz):", (f_hi_ghz - f_lo_ghz) * 1000**2)




compute_parameters(183, 500, 0.05)
