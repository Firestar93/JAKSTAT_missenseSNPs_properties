import matplotlib.pyplot as plt
from matplotlib.patches import Circle
import matplotlib
matplotlib.use('TkAgg')  # Switch to Tkinter backend

venn_diagram_path = "C:\\Users\\hoffmannmd\\OneDrive - National Institutes of Health\\00_PROJECTS\\STAT_JAK_HKs_paper_request\\NucleotideSequence\\EnzymeCuts\\PythonOutput\\BenignAnalysis\\EnzymeOutput\\enzyme_analysis_output\\venn_diagram.png"

# Lists of restriction enzymes
benign_variants = {
    "PspP",
    "FatI",
    "Hin1I",
    "NlaIII",
    "PsuGI",
    "MnlI",
    "RlaI",
    "MspI",
    "EcoO157SI",
    "EcoBLMcrX",
    "FspEI",
    "Hpy188III",
    "SetI",
    "SgeI",
    "SgrTI",
    "YkrI",
    "HpyCH4V",
    "Hpy178III",
    "LpnPI",
    "BmeDI",
    "AspBHI",
    "CviKI_1",
    "CviJI",
    "CviRI",
    "AbaSI"
}

clinical_variants = {
    "Hpy178III", "MnlI", "Lcr047I", "LpnPI", "EcoBLMcrX", "Asi256I", "MalI",
    "Sau3AI", "DpnI", "DpnII", "Bsp143I", "NdeI", "BssMI", "FaiI", "AspBHI",
    "SetI", "SgrTI", "AbaSI", "EcoO157SI", "RlaI", "YkrI", "SgeI", "BmeDI",
    "FspEI", "MspJI"
}

# Calculate subsets
only_benign = benign_variants - clinical_variants
only_clinical = clinical_variants - benign_variants
both = benign_variants & clinical_variants

# Create figure with high DPI
fig, ax = plt.subplots(figsize=(12, 10), dpi=300)  # Increase figure size & DPI
ax.set_xlim(0, 10)
ax.set_ylim(0, 10)

# Draw circles
circle1 = Circle((4, 5), 3.8, edgecolor='#352585', facecolor='none', linewidth=2, label="ClinVar benign/likely benign variants")
circle2 = Circle((6, 5), 3.8, edgecolor='#a22f43', facecolor='none', linewidth=2, label="Disease-associated variants")

ax.add_patch(circle1)
ax.add_patch(circle2)

# Add text labels with higher resolution
ax.text(1.5, 3.5, "\n".join(only_benign), fontsize=12, ha='center', color='#352585', fontweight='bold')
ax.text(8.5, 3.5, "\n".join(only_clinical), fontsize=12, ha='center', color='#a22f43', fontweight='bold')
ax.text(5, 3, "\n".join(both), fontsize=12, ha='center', color='#2b2b2b', fontweight='bold')

# Title and legend
ax.set_title("Venn Diagram of Restriction Enzymes", fontsize=16)
ax.legend(loc="upper right")

# Hide axes
ax.set_xticks([])
ax.set_yticks([])
ax.set_frame_on(False)

# Save as high-resolution image
plt.savefig(venn_diagram_path, dpi=600, bbox_inches='tight')  # Save at 600 DPI

plt.show()
