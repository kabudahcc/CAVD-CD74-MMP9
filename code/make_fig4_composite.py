# -*- coding: utf-8 -*-
"""重画 WGCNA GS-MM 散点（turquoise 模块，高亮 MMP9）并合成新 Fig. 4：
A = 模块-性状热图（R 生成）；B = GS vs MM；C = STRING PPI
"""
import os
import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from PIL import Image, ImageDraw, ImageFont

OUT = r"C:\Users\98757\Desktop\生信分析\论文_选题B"
m = pd.read_csv(os.path.join(OUT, "TableS_WGCNA_module_membership.csv"))
t = m[m["module"] == "turquoise"]

fig, ax = plt.subplots(figsize=(5.2, 4.6), dpi=300)
ax.scatter(t["MM_key"].abs(), t["GS"].abs(), s=10, c="#78909C", alpha=0.6, edgecolors="none")
r = np.corrcoef(t["MM_key"].abs(), t["GS"].abs())[0, 1]
row = m[m["gene"] == "MMP9"].iloc[0]
ax.scatter([abs(row["MM_key"])], [abs(row["GS"])], marker="^", s=120, c="#C0504D",
           edgecolors="black", linewidths=0.8, zorder=5)
ax.annotate("MMP9", (abs(row["MM_key"]), abs(row["GS"])), textcoords="offset points",
            xytext=(10, 6), fontsize=11, fontweight="bold", color="#C0504D")
ax.set_xlabel("Module membership (turquoise)", fontsize=9)
ax.set_ylabel("Gene significance for calcification", fontsize=9)
ax.set_title(f"turquoise module (n = {len(t)}): GS vs MM, r = {r:.2f}", fontsize=10)
ax.tick_params(labelsize=8)
fig.tight_layout()
fig.savefig(os.path.join(OUT, "Fig4_WGCNA_GS_MM.png"), bbox_inches="tight")
plt.close(fig)
print("GS-MM r =", round(r, 3))

# ---- 合成 Fig. 4：A 模块-性状热图 | B GS-MM | C PPI ----
def load_font(sz):
    for p in [r"C:\Windows\Fonts\arialbd.ttf", r"C:\Windows\Fonts\arial.ttf"]:
        if os.path.exists(p):
            return ImageFont.truetype(p, sz)
    return ImageFont.load_default()

def pad_white(img, w, h, label):
    canvas = Image.new("RGB", (w, h), "white")
    img.thumbnail((w - 20, h - 70))
    canvas.paste(img, ((w - img.width) // 2, 60 + (h - 60 - img.height) // 2))
    d = ImageDraw.Draw(canvas)
    d.text((w * 0.02, h * 0.015), label, fill="black", font=load_font(max(30, w // 26)))
    return canvas

imgs = [Image.open(os.path.join(OUT, f)).convert("RGB")
        for f in ("Fig4_WGCNA_module_trait.png", "Fig4_WGCNA_GS_MM.png", "Fig4_PPI_network.png")]
cellw = max(i.width for i in imgs); cellh = max(i.height for i in imgs)
canvas = Image.new("RGB", (cellw * 3, cellh), "white")
for k, (im, lab) in enumerate(zip(imgs, "ABC")):
    canvas.paste(pad_white(im.copy(), cellw, cellh, lab), (k * cellw, 0))
canvas.save(os.path.join(OUT, "Fig4_WGCNA_PPI.png"))
print("saved Fig4_WGCNA_PPI.png", canvas.size)
