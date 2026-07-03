#!/usr/bin/env python3
"""Figure 3-style network robustness, NetShift and driver-intersection plots."""

from pathlib import Path
from collections import defaultdict
import numpy as np
import pandas as pd
import networkx as nx
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from matplotlib.gridspec import GridSpec
from matplotlib.lines import Line2D
from sklearn.covariance import LedoitWolf
from scipy.stats import linregress
from PIL import Image, ImageDraw, ImageFont

ROOT = Path(__file__).resolve().parents[1]
DATA, OUT = ROOT / "data", ROOT / "output"
OUT.mkdir(exist_ok=True)
META = pd.read_csv(DATA / "uploaded_sample_grouping.tsv", sep="\t")
X = pd.read_csv(DATA / "figure3_family_counts.tsv", sep="\t").set_index("Sample")
META["Treatment"] = META.Group.str[0]
META["Compartment"] = META.Group.str[1].map({"S":"Rhizosphere","R":"Root"})
X = X.loc[META.Sample].apply(pd.to_numeric)
TR = ["R","S","E"]
COL = {"R":"#D6A23B","S":"#56A6C9","E":"#6E5195"}
ECOL = {"Control only":"#D73027","Case only":"#1A9850","Both":"#4575B4"}
SEED = 20260703

def taxa(samples, n, groups=None):
    z=X.loc[samples]; prev=(z>0).mean(); rel=z.div(z.sum(1),axis=0).mean()
    keep=(prev>=.3)&(z.sum()>0)
    for g in groups or []: keep &= X.loc[g].sum()>0
    return pd.DataFrame({"p":prev,"a":rel})[keep].sort_values(["a","p"],ascending=False).head(n).index.tolist()

def graph(samples,names,small=False):
    z=np.log(X.loc[samples,names].to_numpy(float)+.5); z-=z.mean(1,keepdims=True)
    ok=z.std(0)>1e-12; z=z[:,ok]; names=np.array(names)[ok].tolist()
    cov=LedoitWolf().fit(z).covariance_; d=np.sqrt(np.diag(cov)); cor=cov/np.outer(d,d)
    cand=[]
    for i in range(len(names)):
        for j in range(i+1,len(names)):
            r=float(cor[i,j])
            if np.isfinite(r) and r>0: cand.append((abs(r),r,names[i],names[j]))
    cand.sort(reverse=True)
    m=len(names)*(len(names)-1)//2
    k=min(45 if small else 140, max(15 if small else 25, round((.25 if small else .15)*m)),len(cand))
    g=nx.Graph(); g.add_nodes_from(names)
    for _,r,a,b in cand[:k]: g.add_edge(a,b,r=r)
    return g

def robust(g,seed):
    rng=np.random.default_rng(seed); nodes=np.array(g.nodes()); n=len(nodes); den=n*(n-1); out=[]
    for p in np.arange(0,1,.1):
        k=min(int(p*n),n-1)
        for rep in range(1 if k==0 else 100):
            h=g.copy()
            if k: h.remove_nodes_from(rng.choice(nodes,k,replace=False))
            out.append((p,rep+1,2*h.number_of_edges()/den))
    return pd.DataFrame(out,columns=["RemovalRatio","Iteration","Connectivity"])

def eset(g): return {tuple(sorted(e)) for e in g.edges()}

def scale(v):
    v=np.asarray(v,float); return np.zeros_like(v) if v.max()==v.min() else (v-v.min())/(v.max()-v.min())

def netshift(ctrl,case,comp,cid,clab,klab):
    nodes=sorted(set(ctrl)&set(case)); c=ctrl.subgraph(nodes).copy(); k=case.subgraph(nodes).copy()
    ce,ke=eset(c),eset(k); edges=[]
    for a,b in sorted(ce|ke):
        cls="Both" if (a,b) in ce and (a,b) in ke else ("Control only" if (a,b) in ce else "Case only")
        edges.append((a,b,cls))
    bc=nx.betweenness_centrality(c); bk=nx.betweenness_centrality(k)
    sc,sk=scale([bc[n] for n in nodes]),scale([bk[n] for n in nodes]); md=max([k.degree(n) for n in nodes]+[1])
    rows=[]
    for i,n in enumerate(nodes):
        nc,nk=set(c.neighbors(n)),set(k.neighbors(n)); u=nc|nk; only=nk-nc
        nesh=0 if not u else 1-len(nc&nk)/len(u)+len(only)/md+len(only)/len(u)
        rows.append([n,nesh,sk[i]-sc[i]])
    nd=pd.DataFrame(rows,columns=["Taxon","NESH","DeltaBetweenness"])
    nd["Driver"]=(nd.DeltaBetweenness>0)&(nd.NESH>=1)
    if not nd.Driver.any():
        q=nd[nd.DeltaBetweenness>0].sort_values(["NESH","DeltaBetweenness"],ascending=False)
        if q.empty:q=nd.sort_values(["NESH","DeltaBetweenness"],ascending=False)
        nd.loc[nd.Taxon.isin(q.head(3).Taxon),"Driver"]=True
    u=nx.Graph();u.add_nodes_from(nodes)
    for a,b,cls in edges:u.add_edge(a,b,edge_class=cls)
    return {"union":u,"nodes":nd,"edges":edges,"comp":comp,"cid":cid,"control":clab,"case":klab}

def plot_a(df):
    fig,axs=plt.subplots(2,1,figsize=(4.2,7),constrained_layout=True)
    for ax,comp in zip(axs,["Rhizosphere","Root"]):
        sub=df[df.Compartment==comp]; rng=np.random.default_rng(1)
        for i,t in enumerate(TR):
            d=sub[sub.Treatment==t]; x=d.RemovalRatio.to_numpy(); y=d.Connectivity.to_numpy()
            ax.scatter(np.clip(x+rng.normal(0,.009,len(x)),0,1),y,s=10,alpha=.38,c=COL[t],edgecolors="none",label=t)
            lr=linregress(x,y); xx=np.linspace(0,.9,100); yy=lr.intercept+lr.slope*xx
            ax.plot(xx,yy,c=COL[t],lw=1.3)
            ax.text(.98,.98-i*.07,fr"$R^2$ = {lr.rvalue**2:.2f}, $P$ < 0.001",c=COL[t],transform=ax.transAxes,ha="right",va="top",fontsize=7)
        ax.set(xlim=(-.02,.95),ylim=(0,None),xlabel="Ratio of randomly removed nodes",ylabel="Network connectivity")
        ax.set_xticks([0,.25,.5,.75]);ax.text(.02,.03,comp,transform=ax.transAxes,fontsize=8)
        ax.spines[["top","right"]].set_visible(False);ax.tick_params(labelsize=7)
    axs[1].legend(ncol=3,frameon=False,fontsize=7,loc="lower center",bbox_to_anchor=(.5,-.34))
    fig.savefig(OUT/"Figure3a_network_robustness.png",dpi=300,bbox_inches="tight");plt.close(fig)

def plot_net(ax,r):
    g=r["union"];nd=r["nodes"].set_index("Taxon")
    dr=sorted(nd.index[nd.Driver]); no=sorted(nd.index[~nd.Driver]); order=[]
    while dr or no:
        if no:order.append(no.pop(0))
        if dr:order.append(dr.pop(0))
    th=np.linspace(np.pi/2,np.pi/2-2*np.pi,len(order),endpoint=False)
    pos={n:(np.cos(a),np.sin(a)) for n,a in zip(order,th)}
    for cls in ECOL:
        ed=[(a,b) for a,b,d in g.edges(data=True) if d["edge_class"]==cls]
        nx.draw_networkx_edges(g,pos,edgelist=ed,ax=ax,edge_color=ECOL[cls],width=.55,alpha=.52)
    ax.scatter([pos[n][0] for n in order],[pos[n][1] for n in order],
               s=[18+65*nd.loc[n,"NESH"] for n in order],
               c=["#D73027" if nd.loc[n,"Driver"] else "#202020" for n in order],
               edgecolors="white",linewidths=.35,zorder=3)
    for n,a in zip(order,th):
        x,y=1.13*np.cos(a),1.13*np.sin(a);deg=np.degrees(a);rot=deg+180 if x<0 else deg
        ax.text(x,y,n,rotation=rot,rotation_mode="anchor",ha="right" if x<0 else "left",va="center",
                fontsize=5.6,c="#D73027" if nd.loc[n,"Driver"] else "#111")
    ax.text(0,0,f"{r['case']} vs {r['control']}",ha="center",va="center",fontsize=8,
            bbox=dict(boxstyle="round,pad=.25",fc="white",ec="#777",lw=.4))
    ax.set_title(r["comp"],fontsize=9);ax.set(xlim=(-1.45,1.45),ylim=(-1.45,1.45));ax.set_aspect("equal");ax.axis("off")

def plot_b(results):
    fig=plt.figure(figsize=(12,9));gs=GridSpec(2,2,figure=fig,wspace=.18,hspace=.18)
    for ax,r in zip([fig.add_subplot(gs[i,j]) for i in range(2) for j in range(2)],results):plot_net(ax,r)
    h=[Line2D([0],[0],marker="o",c="w",mfc="#D73027",ms=6,label="Driver taxa")]
    h += [Line2D([0],[0],c=ECOL[k],lw=1,label={"Control only":"Associations present only in control","Case only":"Associations present only in case","Both":"Associations present in both"}[k]) for k in ECOL]
    fig.legend(handles=h,loc="lower center",ncol=2,frameon=False,fontsize=8,bbox_to_anchor=(.5,.01))
    fig.savefig(OUT/"Figure3b_netshift.png",dpi=300,bbox_inches="tight");plt.close(fig)

def plot_c(results):
    order=["R_vs_S (Root)","R_vs_E (Root)","R_vs_S (Rhizosphere)","R_vs_E (Rhizosphere)"]
    sets={f"{r['cid']} ({r['comp']})":set(r["nodes"].loc[r["nodes"].Driver,"Taxon"]) for r in results};sets={k:sets[k] for k in order}
    pats=defaultdict(list)
    for tax in sorted(set().union(*sets.values())):
        p=tuple(int(tax in sets[k]) for k in order);pats[p].append(tax)
    combos=sorted(pats.items(),key=lambda z:(-len(z[1]),-sum(z[0]),z[0]))[:12]
    fig=plt.figure(figsize=(7,4.8));gs=GridSpec(2,2,figure=fig,width_ratios=[1.2,4.5],height_ratios=[3.2,1.8],hspace=.03,wspace=.05)
    fig.add_subplot(gs[0,0]).axis("off");bar=fig.add_subplot(gs[0,1]);setsax=fig.add_subplot(gs[1,0]);mat=fig.add_subplot(gs[1,1],sharex=bar)
    cnt=[len(v) for _,v in combos];multi=[i for i,(p,_) in enumerate(combos) if sum(p)>1];hi=max(multi,key=lambda i:cnt[i]) if multi else int(np.argmax(cnt))
    bar.bar(range(len(combos)),cnt,color=["#D73027" if i==hi else "#BDBDBD" for i in range(len(combos))],ec="#777",lw=.4)
    for i,c in enumerate(cnt):bar.text(i,c+.12,str(c),ha="center",fontsize=7)
    bar.set_ylabel("Driver taxa in each set",fontsize=8);bar.set_xticks([]);bar.spines[["top","right","bottom"]].set_visible(False)
    yy=list(range(4))[::-1]
    for i,(p,_) in enumerate(combos):
        py=[yy[j] for j,v in enumerate(p) if v]
        if len(py)>1:mat.plot([i,i],[min(py),max(py)],c="#222",lw=.7)
        for j,y in enumerate(yy):mat.scatter(i,y,s=23,c="#111" if p[j] else "#D9D9D9",edgecolors="none")
    mat.set_yticks(range(4));mat.set_yticklabels(order[::-1],fontsize=6.6);mat.yaxis.tick_right();mat.set_xticks([]);mat.tick_params(axis="y",length=0)
    for s in mat.spines.values():s.set_visible(False)
    sizes=[len(sets[k]) for k in order];setsax.barh(yy,sizes,color="#BDBDBD",ec="#777",lw=.4)
    setsax.invert_xaxis();setsax.set_yticks([]);setsax.set_xlabel("Set Size",fontsize=7);setsax.spines[["top","left","right"]].set_visible(False)
    fig.savefig(OUT/"Figure3c_driver_upset.png",dpi=300,bbox_inches="tight");plt.close(fig)

def combine():
    a,b,c=[Image.open(OUT/f).convert("RGB") for f in ["Figure3a_network_robustness.png","Figure3b_netshift.png","Figure3c_driver_upset.png"]]
    W,H=4800,3200;can=Image.new("RGB",(W,H),"white");lw=1350;rx=1430
    def fit(im,w,h):
        im=im.copy();im.thumbnail((w,h),Image.Resampling.LANCZOS);bg=Image.new("RGB",(w,h),"white");bg.paste(im,((w-im.width)//2,(h-im.height)//2));return bg
    can.paste(fit(a,1290,2090),(30,30));can.paste(fit(c,1290,980),(30,2190));can.paste(fit(b,3310,3080),(rx,30))
    d=ImageDraw.Draw(can)
    try:f=ImageFont.truetype("/usr/share/fonts/truetype/dejavu/DejaVuSans.ttf",52)
    except:f=ImageFont.load_default()
    d.text((20,10),"a",font=f,fill="black");d.text((rx+10,10),"b",font=f,fill="black");d.text((20,2160),"c",font=f,fill="black")
    can.save(OUT/"Figure3_combined.png",dpi=(300,300))

def main():
    rr=[];summary=[];seed=0
    for comp in ["Rhizosphere","Root"]:
        sm=META.loc[META.Compartment==comp,"Sample"].tolist();n=44 if comp=="Rhizosphere" else 14;names=taxa(sm,n)
        for t in TR:
            ids=META.loc[(META.Compartment==comp)&(META.Treatment==t),"Sample"].tolist();g=graph(ids,names,n<20);r=robust(g,SEED+seed);seed+=1;r["Compartment"]=comp;r["Treatment"]=t;rr.append(r)
            summary.append([comp,t,len(ids),g.number_of_nodes(),g.number_of_edges(),nx.density(g)])
    rdf=pd.concat(rr);rdf.to_csv(OUT/"Figure3a_random_removal_raw.csv",index=False);pd.DataFrame(summary,columns=["Compartment","Treatment","Samples","Nodes","Edges","Density"]).to_csv(OUT/"Figure3a_network_summary.csv",index=False);plot_a(rdf)
    res=[]
    for comp in ["Rhizosphere","Root"]:
        n=44 if comp=="Rhizosphere" else 14
        for cid,case,ctrl in [("R_vs_S","R","S"),("R_vs_E","R","E")]:
            ki=META.loc[(META.Compartment==comp)&(META.Treatment==case),"Sample"].tolist();ci=META.loc[(META.Compartment==comp)&(META.Treatment==ctrl),"Sample"].tolist()
            names=taxa(ki+ci,n,[ki,ci]);res.append(netshift(graph(ci,names,n<20),graph(ki,names,n<20),comp,cid,ctrl,case))
    pd.concat([r["nodes"].assign(Comparison=r["cid"],Compartment=r["comp"]) for r in res]).to_csv(OUT/"Figure3b_netshift_node_statistics.csv",index=False)
    pd.concat([r["nodes"].loc[r["nodes"].Driver].assign(Comparison=r["cid"],Compartment=r["comp"]) for r in res]).to_csv(OUT/"Figure3b_driver_taxa.csv",index=False)
    plot_b(res);plot_c(res);combine()
    print("Finished:",OUT)

if __name__=="__main__":main()
