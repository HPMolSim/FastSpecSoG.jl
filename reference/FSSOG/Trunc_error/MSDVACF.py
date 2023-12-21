import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
#df = pd.read_csv("tips1.csv")
df = pd.read_csv("NewDataMSD.csv")
print(df)
#ax = sns.violinplot(x="Algorithm", y="Temperature (Kelvin)", data=df)#,inner=True, whis=np.inf)
#plt.xlabel('Algorithm',fontdict={'size':14},labelpad=10)

#ax.set_ylim(290,340)
#ax.set_xticks([])#不显示x轴刻度标注
#a=np.linspace(0,60000,601)
#ax = sns.swarmplot(x="Algorithm", y="Temperature (Kelvin)", data=df.iloc[a,[1,3]], color="c")


###    Figure 2 Cell   ###
#df = pd.read_csv("cell_density.csv")
#print(df)
#plt.figure()
#ax = sns.violinplot(x="Algorithm", y="Density", data=df)#,inner=True, whis=np.inf)
#ax = sns.violinplot(x="Algorithm", y="Cell Kinetic Energy", data=df)#,inner=True, whis=np.inf)
#plt.ylabel('Density [g/cm$^3$]',fontdict={'size':14},labelpad=10)
#plt.tick_params(labelsize=11)#设置刻度字体大小
#ax.set_xlabel('')



#fig,ax_arr = plt.subplots(1,1, figsize=(12,5))

#########################    Figure 3 MSD   ###########################
#plt.figure()

# Set index
df = pd.read_csv("NewDataMSD.csv")
#b= np.arange(48000,60000, 1);
#a1=np.array([1,2,5,10,20,50,100,200,500,1000,2000,5000,10000])
#a2=[x + 12000 for x in a1]
#a3=[x + 12000 for x in a2]
#a4=[x + 12000 for x in a3]
#a5=np.append(a1,a2)
#a6=np.append(a3,a4)
#a5=np.append(a5,a6)
#b= [x + 12000 for x in a4]
#a=np.append(a5,b)

#ax = sns.lineplot(x="t [ps]", y="MSD [nm^2]",hue="Algorithm",style="Algorithm",
#    markers=True, dashes=False, data=df.iloc[a,[0,1,2,3]])#,estimator=None)
ax1=sns.lineplot(x="Term", y="Error",hue="Base",style="Base",
    markers=True, dashes=False, data=df)


ax1.set_yscale('log')
ax1.set_xscale('log')
ax1.set_ylim(1e-12,1e+3)
ax1.set_xlim(0.5,500)
ax1.set_xlabel("Term")
ax1.set_ylabel("Error")
ax1.tick_params(labelsize=11)

# Might need to loop through the list if there are multiple lines on the plot
ax1.lines[0].set_color('tab:blue')
ax1.lines[1].set_color('tab:orange')
ax1.lines[2].set_color('tab:green')
ax1.lines[3].set_color('tab:red')
ax1.lines[4].set_color('tab:purple')


ax1.lines[1].set_linestyle("")
ax1.lines[1].set_marker("s")
ax1.lines[0].set_linestyle("")
ax1.lines[0].set_marker("o")
ax1.lines[2].set_linestyle("")
ax1.lines[2].set_marker("^")
ax1.lines[3].set_linestyle("")
ax1.lines[3].set_marker("D")
ax1.lines[4].set_linestyle("")
ax1.lines[4].set_marker("*")


ax1.lines[5].set_color('tab:blue')
ax1.lines[6].set_color('tab:orange')
ax1.lines[7].set_color('tab:green')
ax1.lines[8].set_color('tab:red')
ax1.lines[9].set_color('tab:purple')

ax1.lines[6].set_linestyle("-.")
ax1.lines[6].set_marker('')
ax1.lines[5].set_linestyle("-.")
ax1.lines[5].set_marker('')
ax1.lines[7].set_linestyle("-.")
ax1.lines[7].set_marker('')
ax1.lines[8].set_linestyle("-.")
ax1.lines[8].set_marker('')
ax1.lines[9].set_linestyle("-.")
ax1.lines[9].set_marker('')

for i in range(0,5):
    ax1.lines[i].set_markersize(12)
    ax1.lines[i].set_markeredgecolor('none')
    ax1.lines[i].set_linewidth(3)


ax1.lines[3].set_markersize(10)
ax1.lines[3].set_markeredgecolor('none')
ax1.lines[3].set_linewidth(3)
ax1.lines[4].set_markersize(15)
ax1.lines[4].set_markeredgecolor('none')
ax1.lines[4].set_linewidth(3)

for i in range(5, 10):
    ax1.lines[i].set_markersize(12)
    ax1.lines[i].set_markeredgecolor('none')
    ax1.lines[i].set_linewidth(1)

ax1.lines[10].set_color('tab:blue')
ax1.lines[11].set_color('tab:orange')
ax1.lines[12].set_color('tab:green')
ax1.lines[13].set_color('tab:red')
ax1.lines[14].set_color('tab:purple')

ax1.lines[11].set_linestyle("--")
ax1.lines[11].set_marker('')
ax1.lines[10].set_linestyle("--")
ax1.lines[10].set_marker('')
ax1.lines[12].set_linestyle("--")
ax1.lines[12].set_marker('')
ax1.lines[13].set_linestyle("--")
ax1.lines[13].set_marker('')
ax1.lines[14].set_linestyle("--")
ax1.lines[14].set_marker('')


#ax1.lines[4].set_linewidth(3)
#ax1.lines[4].set_markersize(16)
#ax1.lines[3].set_markersize(9)

ax1.set_xlabel('$M$',fontdict={'size':16},labelpad=8.5)
ax1.set_ylabel('$\mathscr{R}^{\mathscr{N}}(M)$',fontdict={'size':16},labelpad=5.5)

ax1.legend(loc='best',
          fontsize='large',
          markerscale=1,
          markerfirst=True,
          numpoints=1,
          scatterpoints=1,
          handlelength=2,
          frameon=False)

ax1.legend(labels=["b=2","b=1.630","b=1.488","b=1.321","b=1.218"],title=None,frameon=True, fontsize = '11', title_fontsize = "13",framealpha=1,fancybox=True,shadow=True,borderpad=0.5,
          ncol=3)

#legend=ax.legend()
#legend.texts[-1].set_text("Algorithm")


#plt.show()


#########################    Figure 4 VACF   ###########################
#plt.figure()

# Set index

plt.subplots_adjust(wspace=0.4,top=0.950,bottom=0.15,left=0.170,right=0.950)
plt.show()
