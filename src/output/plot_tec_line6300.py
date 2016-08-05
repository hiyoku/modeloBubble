#!/usr/bin python
import os
import matplotlib as mpl
import numpy as np
import matplotlib.cm as cm
import matplotlib.pyplot as plt


def plot(option):
    x = []
    y = []
    z = []

    if int(option) == 1:
        archive = 'output-tec.dat'
        name = 'Output TEC'
    elif int(option) == 2:
        archive = 'output-line6300.dat'
        name = 'Output Line 6300'
    elif int(option) == 4:
        archive = 'output-den.dat'
        name = 'Output DEN'
    else:
        archive = 'output-line1356.dat'
        name = 'Output Line 1356'

    if int(option) == 4:
        with open(archive, 'r') as e:
            for line in e:
                line = line.split('   ')
                x.append(np.float(float(line[0])))
                y.append(np.float(float(line[1])))
                z.append(np.float(float(line[2])))
    else:
        with open(archive, 'r') as e:
            for line in e:
                line = line.split('   ')
                x.append(np.float(float(line[0])))
                y.append(np.float(float(line[1])))
                z.append(np.float(float(line[2])))

    plt.figure(int(option))

    # plt.subplot(211)
    plt.title(name)
    plt.tricontourf(x, y, z, cmap=cm.autumn, antialiased=True)
    plt.colorbar()
    plt.show()

if __name__ == '__main__':
    while True:
        os.system('clear')
        print('======================== Plot =========================')
        print('[1] - Plotar TEC')
        print('[2] - Plotar Line6300')
        print('[3] - Plotar Line1356')
        print('[4] - Sair')
        print('\n\n')

        try:
            option = input('Digite um numero: ')
        except:
            option = 5

        if int(option) == 4:
            os.system('clear')
            break
        elif 0 < int(option) < 4:
            plot(option)
        else:
            pass
