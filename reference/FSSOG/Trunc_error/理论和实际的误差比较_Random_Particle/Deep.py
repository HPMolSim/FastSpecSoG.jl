import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#df = pd.read_csv("tips1.csv")
df = pd.read_csv("data1.csv")

fig,ax_arr = plt.subplots(1,2, figsize=(11,5))

##############################    Figure A   ###################################
#plt.figure()

# Set index
df = pd.read_csv("data1.csv")
ax1=sns.lineplot(x="x_dot", y="y_dot",hue="base_dot",style="base_dot",
    markers=True, dashes=False, data=df,ax=ax_arr[0])


ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e-11,10)
ax1.set_xlim(0.75,500)
ax1.set_xlabel("x_dot")
ax1.set_ylabel("y_dot")
ax1.tick_params(labelsize=11)

# Might need to loop through the list if there are multiple lines on the plot

ax1.lines[0].set_marker("o")
ax1.lines[0].set_linestyle("")
ax1.lines[0].set_color((0/255, 0/255, 255/255))
#ax1.lines[0].set_markeredgecolor('#6BB952')
ax1.lines[1].set_marker("s")
ax1.lines[1].set_linestyle("")
ax1.lines[1].set_color((0/255, 0/255, 255/255))
#ax1.lines[1].set_markeredgecolor('#EC748B')
ax1.lines[2].set_marker("^")
ax1.lines[2].set_linestyle("")
ax1.lines[2].set_color((0/255, 0/255, 255/255))
#ax1.lines[2].set_markeredgecolor('#C4A751')
ax1.lines[3].set_marker("D")
ax1.lines[3].set_linestyle("")
ax1.lines[3].set_color((0/255, 0/255, 255/255))
#ax1.lines[3].set_markeredgecolor('#6EB1DE')
ax1.lines[4].set_marker("*")
ax1.lines[4].set_linestyle("")
ax1.lines[4].set_color((0/255, 0/255, 255/255))
#ax1.lines[4].set_markeredgecolor('#B67FB3')


for i in range(0,5):
    ax1.lines[i].set_markersize(9)
    #ax1.lines[i].set_markerfacecolor('none')
    #ax1.lines[i].set_markeredgewidth(2)
    ax1.lines[i].set_markeredgecolor('none')

ax1.lines[1].set_markersize(8)
ax1.lines[2].set_markersize(8.5)
ax1.lines[3].set_markersize(7.5)
ax1.lines[4].set_markersize(12.5)


ax1.lines[5].set_linestyle("--")
ax1.lines[5].set_marker("")
ax1.lines[5].set_color((0/255, 0/255, 0/255))
ax1.lines[6].set_linestyle("--")
ax1.lines[6].set_marker("")
ax1.lines[6].set_color((0/255, 0/255, 0/255))
ax1.lines[7].set_linestyle("--")
ax1.lines[7].set_marker("")
ax1.lines[7].set_color((0/255, 0/255, 0/255))
ax1.lines[8].set_linestyle("--")
ax1.lines[8].set_marker("")
ax1.lines[8].set_color((0/255, 0/255, 0/255))
ax1.lines[9].set_linestyle("--")
ax1.lines[9].set_marker("")
ax1.lines[9].set_color((0/255, 0/255, 0/255))

for i in range(5,10):
    ax1.lines[i].set_linewidth(1.5)


ax1.set_xlabel('$M$',fontdict={'size':16},labelpad=8.5)
ax1.set_ylabel('Error of Energy',fontdict={'size':16},labelpad=12.5)

ax1.legend(loc='best',
          fontsize='large',
          markerscale=1,
          markerfirst=True,
          numpoints=1,
          scatterpoints=1,
          handlelength=2,
          frameon=False)

ax1.legend(labels=["$b=2$","$b=1.630$","$b=1.488$","$b=1.321$","$b=1.218$"],title=None,frameon=True, fontsize = '11', title_fontsize = "13",framealpha=1,fancybox=True,shadow=True,borderpad=0.5,loc=3)

#legend=ax.legend()
#legend.texts[-1].set_text("Algorithm")
ax1.minorticks_on()
ax1.text(250,2.74e-1,'A', fontsize=24)

#plt.show()


###########################    Figure B   ##################################

# Set index
df = pd.read_csv("data2.csv")

ax1=sns.lineplot(x="x_dot", y="y_dot",hue="base_dot",style="base_dot",
    markers=True, dashes=False, data=df,ax=ax_arr[1])


ax1.set_yscale('log',subs=[2,3,4,5,6,7,8,9])
ax1.set_xscale('log')
ax1.set_ylim(1e-10,1e+2)
ax1.set_xlim(0.75,500)
ax1.set_xlabel("x_dot")
ax1.set_ylabel("y_dot")
ax1.tick_params(labelsize=11)

# Might need to loop through the list if there are multiple lines on the plot

for i in range(0,5):
    ax1.lines[i].set_markersize(9)
    #ax1.lines[i].set_markerfacecolor('none')
    ax1.lines[i].set_markeredgecolor('none')

ax1.lines[0].set_marker("o")
ax1.lines[0].set_linestyle("")
ax1.lines[0].set_color((0/255, 0/255, 255/255))
#ax1.lines[0].set_markeredgecolor('#6BB952')
ax1.lines[1].set_marker("s")
ax1.lines[1].set_linestyle("")
ax1.lines[1].set_color((0/255, 0/255, 255/255))
#ax1.lines[1].set_markeredgecolor('#EC748B')
ax1.lines[2].set_marker("^")
ax1.lines[2].set_linestyle("")
ax1.lines[2].set_color((0/255, 0/255, 255/255))
#ax1.lines[2].set_markeredgecolor('#C4A751')
ax1.lines[3].set_marker("D")
ax1.lines[3].set_linestyle("")
ax1.lines[3].set_color((0/255, 0/255, 255/255))
#ax1.lines[3].set_markeredgecolor('#6EB1DE')
ax1.lines[4].set_marker("*")
ax1.lines[4].set_linestyle("")
ax1.lines[4].set_color((0/255, 0/255, 255/255))
#ax1.lines[4].set_markeredgecolor('#B67FB3')

ax1.lines[1].set_markersize(8)
ax1.lines[2].set_markersize(8.5)
ax1.lines[3].set_markersize(7.5)
ax1.lines[4].set_markersize(12.5)


ax1.lines[5].set_linestyle("--")
ax1.lines[5].set_marker("")
ax1.lines[5].set_color((0/255, 0/255, 0/255))
ax1.lines[6].set_linestyle("--")
ax1.lines[6].set_marker("")
ax1.lines[6].set_color((0/255, 0/255, 0/255))
ax1.lines[7].set_linestyle("--")
ax1.lines[7].set_marker("")
ax1.lines[7].set_color((0/255, 0/255, 0/255))
ax1.lines[8].set_linestyle("--")
ax1.lines[8].set_marker("")
ax1.lines[8].set_color((0/255, 0/255, 0/255))
ax1.lines[9].set_linestyle("--")
ax1.lines[9].set_marker("")
ax1.lines[9].set_color((0/255, 0/255, 0/255))

for i in range(5,10):
    ax1.lines[i].set_linewidth(1.5)


ax1.set_xlabel('$M$',fontdict={'size':16},labelpad=8.5)
ax1.set_ylabel('Error of Force',fontdict={'size':16},labelpad=12.5)

ax1.legend(loc='best',
          fontsize='large',
          markerscale=1,
          markerfirst=True,
          numpoints=1,
          scatterpoints=1,
          handlelength=2,
          frameon=False)

ax1.legend(labels=["$b=2$","$b=1.630$","$b=1.488$","$b=1.321$","$b=1.218$"],title=None,frameon=True, fontsize = '11', title_fontsize = "13",framealpha=1,fancybox=True,shadow=True,borderpad=0.5)

#legend=ax.legend()
#legend.texts[-1].set_text("Algorithm")

ax1.text(250,2,'B', fontsize=24)


ax1.minorticks_on()

plt.subplots_adjust(wspace=0.3,hspace=0.271,top=0.980,bottom=0.15,left=0.100,right=0.983)
plt.show()
