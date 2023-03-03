import matplotlib.pyplot as plt

x = [2, 4, 8, 16, 32, 64, 128, 256] # x-axis
y1_p2p = [0.000303, 0.001116, 0.012199, 0.038164, 0.389912, 0.773272, 0.709906, 0.664659] # y-axis1
y2 = [0.033976, 0.000133, 0.005236, 0.007273, 0.273762, 0.304891, 0.176335, 0.097598] # y-axis2

plt.plot(x, y1_p2p, 'r', label='MPI_P2P_Reduce')
plt.plot(x, y2, 'b', label='MPI_Reduce')
plt.legend(loc='upper left')
plt.xlabel('Number of Ranks')
plt.ylabel('Runtime (s)')
plt.title('MPI_Reduce vs MPI_P2P_Reduce')
plt.show()