import numpy as np
import matplotlib.pyplot as plt


INTEGRATION_TIME = 0.05
N_TIMESTEPS = 100
CENTRE_FREQ = 158
FREQ_RES = 0.01
N_CHANNELS = int(1.28 / FREQ_RES)
K = 4.15 #8808e3

def delay(DM, f1, f2):
    return K * DM * (f1**(-2) - f2**(-2))

# noise = np.random.normal(0, 0.5, (N_CHANNELS, N_TIMESTEPS))

# times = [INTEGRATION_TIME * i for i in range(N_TIMESTEPS)]
# freqs = np.linspace(CENTRE_FREQ - 1.28 / 2, CENTRE_FREQ + 1.28 / 2, N_CHANNELS)
# F = [(0.07 + i * 0.03) for i in range(5)] # GHz 
# v = [delay(500, F[i - 1], F[i]) / 1000 for i in range(1, len(F))]


# plt.plot(v, F[1:])
# plt.xlabel("Delay (s)")
# plt.ylabel("Frequency (GHz)")
# plt.show()
# T = 0.30
# for i in range(N_CHANNELS - 1):
#     print(freqs[i])
#     D = delay(500, freqs[i + 1], freqs[i])
#     print(D)
#     T_adj = T + D
#     #print(T_adj)
#     noise[i, int(T_adj/INTEGRATION_TIME)] = 3
#     #break
# plt.imshow(noise)
# plt.show()