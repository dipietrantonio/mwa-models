from basics import dispersive_delay_ms

def compute_parameters(centre_fq_mhz, dm, tres):
    # TODO: compute for the highest frequency channel in the band?
    cf_ghz = centre_fq_mhz / 1000
    bwghz = 0.03072
    f_lo_ghz = cf_ghz - bwghz / 2
    
    delta_ghz = 0.00001 # 10 kHz
    f_hi_ghz = f_lo_ghz + delta_ghz
    
    delay = dispersive_delay_ms(dm, f_lo_ghz, f_hi_ghz) / 1000
    while delay < tres:
        f_hi_ghz += delta_ghz
        delay = dispersive_delay_ms(dm, f_lo_ghz, f_hi_ghz) / 1000
    
    print("Final bandwidth required (kHz):", (f_hi_ghz - f_lo_ghz) * 1000**2)




compute_parameters(150, 600, 0.05)
