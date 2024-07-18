import numpy as np
from pylab import*
import matplotlib.pyplot as plt  
from commpy.sequences import pnsequence
from scipy.signal import butter, lfilter


def butter_lowpass(cutoff, fs, order=5):
    return butter(order, cutoff, fs=fs, btype='low', analog=False)

def butter_lowpass_filter(data, cutoff, fs, order=5):
    b, a = butter_lowpass(cutoff, fs, order=order)
    y = lfilter(b, a, data)
    return y

def bitfield(n):
    return [int(digit) for digit in bin(n)[2:]] 

def modulate(signal, T_chip, f_carrier, accuracy):
    t = r_[0:T_chip:T_chip/accuracy]
    carrier = sin(2 * pi * f_carrier * t)
    bpskarray=[]
    for i in range(0,len(signal)):
        if signal[i]>=0:
            bpskarray.append(carrier)
        else:
            bpskarray.append((-1*carrier))
    return concatenate(bpskarray)

def demodulate(signal, T_chip, f_carrier, accuracy):
    t = r_[0:T_chip:T_chip/accuracy]
    carrier = sin(2 * pi * f_carrier * t)
    carrier = np.resize(carrier, len(signal))
    return signal * carrier

def detect(despreaded_signal, codelegnth):
    l = int(len(despreaded_signal)/codelegnth)
    signal = despreaded_signal
    recovered = []
    for i in range(l):
        x = 1.0*sum(signal[i*codelegnth:(i+1)*codelegnth])/codelegnth
        recovered = np.append(recovered, np.abs(x))
    recovered = np.repeat(recovered, codelegnth)
    return recovered

def circular_shift(signal, shift, left):
    res = []
    if left:
        res = np.resize(signal, len(signal) + shift)[shift:]
    else:
        temp = np.append(signal, signal)
        res = np.resize(temp[len(signal)-shift:], len(signal))
    return res


def pn_code(gd, shift):
    state = ''.join(map(str,gd))
    code_len = (2 ** len(gd) - 1) 
    res = []
    if shift == 0:
        shift = code_len
    for i in range(shift):
        res = pnsequence(len(gd), state, g30, code_len)
        m = list(res[1:len(gd) + 1]).__reversed__()
        state = ''.join(map(str,m))
    return res


##########################################  Parameters ########################################

# creating G(D) for PN-code
g30 = np.array([1, 0, 0, 0, 0, 0, 1])
codelength = len(g30)
pn_code_len = (2 ** len(g30) - 1) 
code = pn_code(g30, 0)

total_time = 4
f_data = 5
f_chip = 50
f_carrier = 50
accuracy = 5
threshold = 0.7
transmitted_signal_delay = 0 # samples
delta_T = 1 # samples
lables = True
filter_order = 6
filter_fs = f_carrier    # sample rate, Hz
filter_cutoff = 2 * f_data 

proccessing_gain = int(f_chip / f_data)
data_points = f_data * total_time
all_chips = int(f_chip * total_time)
all_points = int(all_chips * accuracy)


##########################################   Transmitter ########################################

# Create data 
data = 0x9f27
data = np.array(bitfield(data))
data = np.resize(data, data_points)
data = np.repeat(data, proccessing_gain)
data = np.resize(data, all_chips)

# Create code
code = np.repeat(code, accuracy)
code = np.resize(code, all_chips)

# spreading data
spread = (data*2-1)*(code*2-1)
despread = (spread)*(code*2-1)

# mixed with code signal
transmitted_signal = modulate(spread, 1/f_carrier, f_carrier, accuracy)

# Creating Delay
transmitted_signal = circular_shift(transmitted_signal, transmitted_signal_delay, True)

##########################################  Transmitter END ########################################

# despreading with the same code to see if it works

data = np.repeat(data*2-1, accuracy)
code = np.repeat(code*2-1, accuracy)
code = circular_shift(code, transmitted_signal_delay, True)

spread = np.repeat(spread, accuracy)
despread = np.repeat(despread, accuracy)

# despreading signal
despreaded = code * transmitted_signal

# demodulated_signal
demodulated = demodulate(despreaded, 1/f_carrier, f_carrier, accuracy)

# Low pass filter
y = butter_lowpass_filter(demodulated, filter_cutoff, filter_fs, filter_order)

detected = detect(y, int(all_points / data_points)) * 2

detected_energy = np.mean(detected)

##########################################   Reciever  ########################################

max_detected_energy = 0
max_detected_energy_code = []
max_detected_energy_despreaded = []
max_detected_energy_demodulated = []
max_detected_energy_y = []
max_detected_energy_detected = []

# start of acquisition loop :
for i in range(pn_code_len):
    code_2 = pn_code(g30, i + 1 % pn_code_len) * 2 - 1
    code_2 = np.repeat(code_2, accuracy)
    code_2 = np.resize(code_2, all_chips)
    code_2 = np.repeat(code_2, accuracy)
    code_2 = circular_shift(code_2, transmitted_signal_delay + delta_T, True)

    # despreading signal
    despreaded_2 = code_2 * transmitted_signal

    # demodulated_signal
    demodulated_2 = demodulate(despreaded_2, 1/f_carrier, f_carrier, accuracy)

    # Low pass filter
    y_2 = butter_lowpass_filter(demodulated_2, filter_cutoff, filter_fs, filter_order)

    detected_2 = detect(y_2, int(all_points / data_points)) * 2

    detected_energy_2 = np.mean(detected_2)
    print(detected_energy_2)
    if detected_energy_2 > max_detected_energy:
        max_detected_energy = detected_energy_2
        max_detected_energy_code = code_2
        max_detected_energy_despreaded = despreaded_2
        max_detected_energy_demodulated = demodulated_2
        max_detected_energy_y = y_2
        max_detected_energy_detected = detected_2

    if detected_energy_2 > threshold:
        print("\nhit after trying :", i+1,"houses")
        break
    elif i == pn_code_len - 1:
        print("\nnot detected anything after trying :", i+1,"houses")


code_2 = max_detected_energy_code
despreaded_2 = max_detected_energy_despreaded
demodulated_2 = max_detected_energy_demodulated
y_2 = max_detected_energy_y
detected_2 = max_detected_energy_detected
detected_energy_2 = max_detected_energy

##########################################   Reciever END ########################################

print("max detected energy by acquisition : ", max_detected_energy , "\n")
print("detected energy with correct code and delay: "  ,detected_energy)
print("detected energy with acquisition : " ,detected_energy_2, "\n")

g, (ax0, ax1,ax2, ax3, ax4, ax5, ax6, ax7, ax8) = plt.subplots(9, sharex=True, sharey=True)

ax0.step(range(all_points), data, color='red')
ax0.axhline(color="purple", linestyle="dashed")

ax1.step(range(all_points), code, color='red')
ax1.step(range(all_points), code_2, color='purple')
ax1.axhline(color="gray", linestyle="dashed")

ax2.step(range(all_points), spread, color="purple")

ax3.step(range(all_points), despread, color="red")
ax3.axhline(color="gray", linestyle="dashed")

ax4.plot(range(all_points),transmitted_signal, color="purple")
ax4.axhline(color="gray", linestyle="dashed")

ax5.plot(range(all_points),despreaded, color="red")
ax5.plot(range(all_points),despreaded_2, color="purple")
ax5.axhline(color="gray", linestyle="dashed")

ax6.plot(range(all_points),demodulated, color="red")
ax6.plot(range(all_points),demodulated_2, color="purple")
ax6.axhline(color="gray", linestyle="dashed")

ax7.plot(range(len(y)),y, color="red")
ax7.plot(range(len(y_2)),y_2, color="purple")
ax7.axhline(color="gray", linestyle="dashed")

ax8.plot(range(len(detected)),detected, color="red")
ax8.plot(range(len(detected_2)),detected_2, color="purple")
ax8.axhline(color="gray", linestyle="dashed")

if lables:
    plt.text(0, 0.5, "detected", size=10, rotation=0.,
         ha="center", va="center",
         bbox=dict(boxstyle="round",
                   ec=(1., 0.5, 0.5),
                   fc=(1., 0.8, 0.8),
                   )
         )
    plt.text(0, 3, "LPF", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 5.5, "demodulated", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 8, "despreded", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 10.5, "transmitted", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 13, "despreaded", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 15.2, "spreaded", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 17.5, "codes", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )
    plt.text(0, 19.7, "data", size=10, rotation=0.,
            ha="center", va="center",
            bbox=dict(boxstyle="round",
                    ec=(1., 0.5, 0.5),
                    fc=(1., 0.8, 0.8),
                    )
            )

g.subplots_adjust(hspace=0.1)
plt.show()