import matplotlib.pyplot as plt
import seaborn
plt.rcParams['figure.figsize'] = (14.0,8.0) # change figure size

def histogram_orbital(X,label,path):
    plt.figure()
    plt.hist(X)
    plt.xlabel(label)
    plt.savefig(path)
