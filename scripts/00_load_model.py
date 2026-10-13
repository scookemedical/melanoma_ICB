"""Load Human-GEM SBML once and cache as pickle for fast reuse."""
import pickle, sys, time
import cobra
src, dst = sys.argv[1], sys.argv[2]
t0 = time.time()
m = cobra.io.read_sbml_model(src)
print(f"loaded in {time.time()-t0:.0f}s", flush=True)
print(f"reactions={len(m.reactions)} metabolites={len(m.metabolites)} genes={len(m.genes)}")
print("objective:", m.objective.expression)
with open(dst, "wb") as fh:
    pickle.dump(m, fh)
print("cached ->", dst)
