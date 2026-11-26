import decoupler as dc, pandas as pd, traceback
try:
    df = dc.op.collectri(organism="human")
    print(type(df), getattr(df, "shape", None)); print(df.head()); print(df.columns.tolist())
    df.to_csv("out/collectri.csv", index=False)
except Exception:
    traceback.print_exc()
