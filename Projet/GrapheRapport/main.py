import numpy as np
import matplotlib.pyplot as plt
from scipy.optimize import curve_fit


def read_file(filename):
    with open(filename, 'r') as file:
        lst = file.readlines()
        abs = np.zeros(len(lst))
        ord = np.zeros(len(lst))
        for i, v in enumerate(lst):
            _, a, _, o = v.split(":")
            abs[i] = float(a)
            ord[i] = float(o)
        return abs, ord


def plot_couple_angle():
    plt.figure()
    theta, couple = read_file("couple/switch.txt")
    plt.plot(theta * 180 / np.pi, couple, color="black", linewidth=4.0)
    mean = sum(couple)/len(couple)
    plt.plot(theta*180/np.pi, np.array([mean]*len(theta)), linestyle="--", label="couple moyen", color="orange")
    for file, color, l in [("bobine2.txt", 'g', 'A'), ("bobine1.txt", 'r', 'B'), ("bobine3.txt", 'b', 'C')]:
        theta, couple = read_file(file)
        plt.plot(theta * 180 / np.pi, couple, color=color, label=l)
    plt.xlabel(r"$\theta [°]$")
    plt.ylabel(r"$couple [N.m]$")
    plt.grid()
    plt.legend()
    plt.show()


def read_file_speed(filename):
    with open(filename, 'r') as file:
        lst = file.readlines()
        abs = np.zeros(len(lst))
        ord = np.zeros(len(lst))
        add = 0
        for i, v in enumerate(lst):
            print(i)
            if i%100 == 0 and i!=0:
                add+=1
            _, a, _, o = v.split(":")
            abs[i] = float(a)+add
            ord[i] = float(o)
        return abs, ord


def plot_vitesse_temps():
    plt.figure()
    lst = [('vitesseang400.txt', 'r'), ('vitesseang838.txt', 'g'), ('vitesseang1667.txt', 'b'),
           ('vitesseang4424.txt', 'm'), ('vitesseang14608.txt', 'orange')]
    for file, color in lst:
        t, v = read_file_speed("vitesse/"+file)
        plt.plot(t, v, color=color, label=file.split('g')[1].split('.')[0])
    plt.xlabel(r"$t [s]$")
    plt.ylabel(r"$\omega [rad/s]$")
    plt.grid()
    plt.legend()
    plt.show()


def plot_temps_maillage():
    plt.figure()
    lst = [('time400.txt', 'r'), ('time838.txt', 'g'), ('time1667.txt', 'b'),
           ('time4424.txt', 'm'), ('time14608.txt', 'orange')]
    for file, color in lst:
        t, v = read_file("timeY/"+file)
        plt.semilogy(t, v*1000, color=color, label=file.split('e')[1].split('.')[0])
    plt.xlabel("itérations [-]")
    plt.ylabel("temps [ms]")
    plt.grid(True, which="both", ls="-")
    plt.legend()
    plt.show()


def plot_renum():
    n = np.array([400, 838, 1667, 4424, 14608])
    meanX = [0, 0, 0, 0, 0]
    meanY = [0, 0, 0, 0, 0]
    meanNO = [0, 0, 0, 0]
    for i, f in enumerate(["time400X.txt", "time838X.txt", "time1667X.txt", "time4424X.txt", "time14608X.txt"]):
        t, v = read_file("timeX/" + f)
        meanX[i] = np.mean(v)

    for i, f in enumerate(["time400.txt", "time838.txt", "time1667.txt", "time4424.txt", "time14608.txt"]):
        t, v = read_file("timeY/" + f)
        meanY[i] = np.mean(v)

    for i, f in enumerate(["time400NO.txt", "time838NO.txt", "time1667NO.txt", "time4424NO.txt"]):
        t, v = read_file("timeNO/" + f)
        meanNO[i] = np.mean(v)

    plt.figure()
    plt.plot(n, meanX, color='r', marker="s", label='FEM_XNUM')
    plt.plot(n, meanY, color='b', marker="s", label='FEM_YNUM')
    plt.plot(n[0:4], meanNO, color='g', marker="s", label='FEM_NO')
    plt.xlabel("taille maillage[-]")
    plt.ylabel("temps [s]")
    plt.xticks(n, n)
    plt.grid()
    plt.legend()
    plt.show()


if __name__ == '__main__':
    ref = -0.045267
    n = np.array([400, 838, 1667, 4424, 14608])
    y = np.array([-0.011850, -0.027798, -0.038487, -0.044461, -0.045267])
    plt.figure()
    plt.plot(n, abs(y-ref)/abs(ref), color='r', marker='s')
    plt.xlabel("taille maillage[-]")
    plt.ylabel("|erreur relative| [-]")
    plt.xticks(n, n)
    plt.tight_layout()
    plt.grid()
    plt.show()
    plot_temps_maillage()