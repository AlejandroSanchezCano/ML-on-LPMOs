'''
import matplotlib.pyplot as plt
import numpy as np

def onpick(event):
    print(event)

fig, ax = plt.subplots()
fig.canvas.mpl_connect('pick_event', onpick)
ax.scatter(np.random.random(10), np.random.random(10), picker = True)
plt.show()
'''

%matplotlib inline

import ipywidgets as widgets
from IPython.display import display
import seaborn as sns
import matplotlib.pyplot as plt
sns.set_style("darkgrid")

dropdown = widgets.Dropdown(
    options=['2009', '2010', '2011', "2012", "2013", "2014", "2015", "2016", "2017", "2018","2019"],
    value='2009',
    description='Jahr:',
    )

def Dropdown_Menu(b):

    fig, ax = plt.subplots()
    fig.dpi = 500

    x = [2009, 2010, 2011, 2012, 2013, 2014, 2015, 2016, 2017, 2018, 2019]
    y = [130, 137, 104, 147, 401, 274, 234, 770, 857, 746, 704]

    ax.plot(x, y, label = "Flugstunden pro Jahr", marker = ".")


    ax.legend()

    ax.set_title("Flugstunden")
    ax.set_xlabel("Jahr")
    ax.set_ylabel("Flugstunden")
    ax.set_facecolor((0.9,0.9,0.9))
    display(dropdown)



plt.show()

dropdown.observe(Dropdown_Menu, names="value")    
display(dropdown)