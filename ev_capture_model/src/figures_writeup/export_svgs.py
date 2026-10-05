"""Re-export figures as SVG with genuinely editable text (svg.fonttype='none')."""
import os, sys, importlib, matplotlib
matplotlib.use('Agg')
matplotlib.rcParams['svg.fonttype'] = 'none'   # keep <text>, do not convert to <path>
HERE=os.path.dirname(os.path.abspath(__file__)); sys.path.insert(0,HERE)
FIG=os.path.abspath(os.path.join(HERE,'..','..','..','Poster, Paper, and Write Up','Figures'))
mods = sys.argv[1:] or [f[:-3] for f in sorted(os.listdir(HERE)) if f.startswith('fig') and f.endswith('.py')]
for name in mods:
    m=importlib.import_module(name)
    out=os.path.join(FIG, name+'.svg')
    # each figure module exposes exactly one plotting entry point; find it
    fn=None
    for cand in ('build','plot_reordering_example','plot_skew_scatter','plot_comparison'):
        if hasattr(m,cand): fn=getattr(m,cand); break
    if fn is None: print('  SKIP (no entry point):', name); continue
    if fn.__name__=='plot_comparison': fn(m.build_comparison(), out)
    else: fn(out)
    n_text=open(out,encoding='utf-8').read().count('<text')
    print('  %-50s <text> elements: %d' % (name+'.svg', n_text))
