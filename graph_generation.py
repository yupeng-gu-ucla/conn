# graph generation

import matplotlib.pyplot as plt
from copy import deepcopy
import random, math, sys
import numpy as np

N = 300     # num of nodes
L = 5       # lattice length
v = 0.03    # velocity
eta = 0     # noise (theta)
MAX_ITER = 200


def aver_vec(thetas):
    sum_t = np.sum(thetas, axis=0) / len(thetas)
    return np.linalg.norm(sum_t)


def build_hard_graph(x, b, TH):
    # build graph using hard threshold TH
    N = len(x)
    #THRESHOLD = math.log(TH/(1.0-TH))
    #THRESHOLD = TH*TH

    '''
    res = np.dot(x, x.T) + b
    for i in range(N):
        res[i,i] = 0.0
    '''

    #'''
    res = np.zeros((N,N))
    for i in range(N):
        res[i,i] = -990.0
        for j in range(i+1, N):
            #d = x[i][0] * x[j][0] + x[i][1] * x[j][1]
            dis = ( (x[i][0] - x[j][0]) * (x[i][0] - x[j][0]) + (x[i][1] - x[j][1]) * (x[i][1] - x[j][1]) )
            d = -dis / b[i] * b[j]
            res[i,j] = res[j,i] = d
    #'''

    res = (res > TH)
    #res = (res > THRESHOLD)
    den = 1.0 * np.sum(res) / (N*(N-1))
    return den, res

def propagate(x, b, v, eta, TH, thetas = []):
    # calculate G^{t}
    den, graph = build_hard_graph(x, b, TH)
    graph = graph + 0.0

    N = len(x)

    if len(thetas) != 0:
        # update x^{t+1} (based on theta^{t})
        x = x + v * thetas

    # calculate new theta^{t+1}
    for n in range(N):
        graph[n] *= np.linalg.norm(x[n])        # weighted
    new_thetas = graph.dot(thetas)
    row_sum = np.sum(graph, axis=1)
    tmp_thetas = np.random.rand(N,1) * 2 * math.pi
    for n in range(N):
        if row_sum[n] != 0:
            tmp_thetas[n,0] = math.atan2(new_thetas[n,1], new_thetas[n,0])
        tmp_thetas[n,0] += (random.random()-0.5) * eta   # noise
    new_thetas = np.concatenate((np.cos(tmp_thetas), np.sin(tmp_thetas)), axis=1)   # N*2, each row is unit vector

    return x, new_thetas, den


def plot_dir(x1, b1, thetas, filename):
    N = len(x1)
    fig = plt.figure(figsize = (8,6))
    ax = plt.axes()
    for n in range(N):
        _x1, _y1 = x1[n][0], x1[n][1]
        _x2, _y2 = _x1 + thetas[n,0], _y1 + thetas[n,1]
        if b1[n] > 1.9:
            plt.plot([_x1, _x2], [_y1, _y2], 'r')
            #plt.scatter(_x2, _y2, edgecolor = 'none', s = 15, c = 'r')
            ax.arrow(_x1, _y1, thetas[n,0],  thetas[n,1], head_width=0.2, head_length=0.12, fc='r', ec='r')
        else:
            plt.plot([_x1, _x2], [_y1, _y2], 'b')
            #plt.scatter(_x2, _y2, edgecolor = 'none', s = 15, c = 'b')
            ax.arrow(_x1, _y1, thetas[n,0], thetas[n,1], head_width=0.2, head_length=0.12, fc='b', ec='b')
    plt.xlim([-7,7]); plt.ylim([-7,7])
    #plt.show()
    plt.savefig(filename, filetype = 'pdf')
    plt.close()


def main(TH, eta):
    #b = np.zeros((N,1))                    # zero bias
    #b = (np.random.rand(N,1) - 0.5) * 4.0  # bias for CONN_dot
    b = (np.random.rand(N,1)) * 1.0 + 1     # bias for CONN
    x = (np.random.rand(N,2) - 0.5) * L
    tmp_thetas = np.random.rand(N,1) * 2 * math.pi
    thetas = np.concatenate((np.cos(tmp_thetas), np.sin(tmp_thetas)), axis=1)   # N*2, each row is unit vector

    aver_v = 0.0; aver_den = 0.0

    num_iter = 0
    out_filename = './gravity/sim_' + str(L) + '_' + str(TH) + '_' + str(eta) + '.pdf'
    while True:
        if num_iter > MAX_ITER:
            break

        x, thetas, den = propagate(x, b, v, eta, TH, thetas)

        if num_iter == 0: print 'Iter, Aver_v, density'
        if (num_iter+1) % 10 == 0: print num_iter, aver_vec(thetas), den
        if (num_iter) % 20 == 0: plot_dir(x, b, thetas, out_filename)

        if num_iter > MAX_ITER-10 and num_iter <= MAX_ITER:
            aver_v += aver_vec(thetas)
            aver_den += den
        num_iter += 1
    #print 'Average v = %f' % (aver_v/10)
    #print 'Average density = %f' % (aver_den/10)
    #print str(L) + ' ' + str(aver_v/10) + ' ' + str(aver_den/10)

    '''
    with open('./gravity/log/aver_vec.log', 'a') as fp:
        newline = str(TH) + '\t' + str(L) + '\t' + str(eta) + '\t' + str(aver_v/10) + '\t' + str(aver_den/10) + '\n'
        fp.write(newline)
    '''


# CONN_dot
def inner_product():
    for _t in [k/4.0 for k in range(1,11)]:
        TH = _t
        #THRESHOLD = math.log(TH/(1.0-TH))
        for e in range(24):
            eta = 1.0*e/4
            print 'TH =', TH, ' eta =', eta
            main(TH, eta)


# CONN
def gravity():
    for _t in [k/4.0 for k in range(1,11)]:
        TH = -_t
        for eta in [l/2.0 for l in range(9)]:
            print 'TH =', TH, ' eta =', eta
            main(TH, eta)


if __name__ == '__main__':
    np.random.seed(42)

    #inner_product()
    gravity()

