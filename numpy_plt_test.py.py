import numpy as np
import matplotlib.pyplot as plt

def sigmoid(x):
    """
    culculate sigmoid 
    """
    return 1 / (1 + np.exp(-x))

def main():
    # a series of sigmoid
    x_values = np.linspace(-10, 10, num=100)
    y_values = sigmoid(x_values)

    
    plt.plot(x_values, y_values)
    plt.xlabel('X Value')
    plt.ylabel('Y Value')
    plt.title('Sigmoid Function')
    plt.show()

if __name__ == '__main__':
    main()