import pandas as pd
import matplotlib.pyplot as plt
import seaborn as sns
import os
import re
from sklearn.preprocessing import MinMaxScaler

# Add this function to sanitize filenames
def sanitize_filename(name):
    # Replace characters that are problematic for filenames
    for char in ['/', '\\', ':', '*', '?', '"', '<', '>', '|', '(', ')']:
        name = name.replace(char, '_')
    # Remove any multiple consecutive underscores
    name = re.sub('_+', '_', name)
    return name.strip()

# Define a function to generate heatmap for any protein
def generate_protein_heatmap(data, protein_name, protein_length, label_y_offset=None):
    # Subset the data for the protein (disable regex to avoid pattern warnings)
    protein_data = data[data["Protein names"].str.contains(protein_name, na=False, case=False, regex=False)]

    peptide_length = 56
    overlap = 28

    expected_peptides = []
    start = 1
    while start <= protein_length:
        end = min(start + peptide_length - 1, protein_length)
        expected_peptides.append((start, end))
        start += (peptide_length - overlap)

    expected_peptides_df = pd.DataFrame(expected_peptides, columns=["start", "end"])
    filtered_data = pd.merge(expected_peptides_df, protein_data, on=["start", "end"], how="inner")

    if filtered_data.empty:
        print(f"⚠️ No matching peptides found for {protein_name}. Skipping.")
        return

    averaged_data = filtered_data.groupby(["start", "end"]).mean(numeric_only=True).reset_index()
    averaged_data["unique_id"] = averaged_data["start"].astype(str) + "_" + averaged_data["end"].astype(str)
    averaged_data = averaged_data.sort_values(by=["start", "end"])
    averaged_data.set_index("unique_id", inplace=True)

    sample_columns = averaged_data.select_dtypes(include=["float64", "int64"]).columns
    heatmap_data = averaged_data[sample_columns].fillna(0)

    print(f"🔍 {protein_name}: {len(heatmap_data)} peptides, {heatmap_data.shape[1]} samples")
    positives = heatmap_data.sum(axis=1)
    if not positives.empty:
        print(f"Most reactive peptide: {positives.idxmax()} ({positives.max()} positives)")
    total_positives = heatmap_data.values.sum()
    total_cells = heatmap_data.size
    if total_cells > 0:
        print(f"Positivity rate: {total_positives}/{total_cells} ({total_positives / total_cells:.2%})")

    # Apply MinMaxScaler to normalize data to 0-1 range (added from original code)
    scaler = MinMaxScaler(feature_range=(0, 1))
    scaled_heatmap_data = pd.DataFrame(
        scaler.fit_transform(heatmap_data),
        index=heatmap_data.index,
        columns=heatmap_data.columns
    )

    sample_groups = {
        "Acute(n = 97)": ["1527v1", "1531v1", "1536v1", "1537v1", "1538v1", "1539v1", "1541v1",
            "1542v1", "1543v1", "1545v1", "1547v1", "1548v1", "1550v1", "1555v1",
            "1556v1", "1558v1", "1559v1", "1563v1", "1564v1", "1567v1", "1568v1",
            "1569v1", "1570v1", "1572v1", "1573v1", "1576v1", "1577v1", "1578v1",
            "1580v1", "1582v1", "1583v1", "1584v1", "1585v1", "1586v1", "1587v1",
            "1588v1", "1589v1", "1590v1", "1591v1", "1594v1", "1598v1", "1600v1",
            "1602v1", "1604v1", "1605v1", "1607v1", "1611v1", "1612v1", "1614v1",
            "1617v1", "1618v1", "1620v1", "1621v1", "1622v1", "1623v1", "1625v1",
            "1626v1", "1627v1", "1630v1", "1633v1", "1634v1", "1635v1", "1636v1",
            "1637v1", "1639v1", "1641v1", "1645v1", "1646v1", "1648v1", "1649v1",
            "1652v1", "1653v1", "1654v1", "1657v1", "1658v1", "1659v1", "1663v1",
            "1664v1", "1665v1", "1666v1", "1667v1", "1670v1", "1671v1", "1674v1",
            "1675v1", "1676v1", "1679v1", "1681v1", "1684v1", "1685v1"],
        "6W(n = 67)": ["1537v5", "1564v5", "1569v5", "1573v5", "1576v5", "1587v5", "1594v5",
            "1611v5", "1618v5", "1623v5", "1625v5", "1652v5", "1653v5", "1666v5",
            "1676v5", "1527v6", "1531v6", "1542v6", "1543v6", "1545v6", "1548v6",
            "1550v6", "1558v6", "1559v6", "1563v6", "1567v6", "1568v6", "1570v6",
            "1575v6", "1577v6", "1582v6", "1585v6", "1588v6", "1590v6", "1600v6",
            "1607v6", "1614v6", "1617v6", "1620v6", "1622v6", "1624v6", "1626v6",
            "1630v6", "1632v6", "1633v6", "1634v6", "1635v6", "1636v6", "1637v6",
            "1639v6", "1641v6", "1645v6", "1646v6", "1648v6", "1649v6", "1651v6",
            "1655v6", "1657v6", "1658v6", "1664v6", "1665v6", "1667v6", "1670v6",
            "1674v6", "1675v6", "1679v6", "1681v6", "ES258v6", "ES308v6",
            "ES405v6", "ES413v6"],
        "6M(n = 30)": ["1531v7", "1538v7", "1563v7", "1564v7", "1568v7", "1577v7", "1578v7",
            "1582v7", "1589v7", "1602v7", "1612v7", "1621v7", "1624v7", "1627v7",
            "1648v7", "1649v7", "1655v7", "1664v7", "1665v7", "1666v7", "1667v7",
            "1670v7", "1671v7", "1674v7", "1675v7", "1676v7", "1679v7", "1681v7",
            "1684v7", "1685v7", "ES308v7", "ES413v7", "ES421v7"],
        "1Y(n = 67)": ["1527v8", "1536v8", "1537v8", "1539v8", "1541v8", "1542v8", "1543v8",
            "1545v8", "1547v8", "1548v8", "1550v8", "1555v8", "1556v8", "1558v8",
            "1559v8", "1567v8", "1569v8", "1570v8", "1572v8", "1573v8", "1575v8",
            "1576v8", "1580v8", "1583v8", "1584v8", "1585v8", "1586v8", "1587v8",
            "1588v8", "1590v8", "1591v8", "1593v8", "1594v8", "1598v8", "1600v8",
            "1603v8", "1604v8", "1605v8", "1607v8", "1611v8", "1614v8", "1617v8",
            "1618v8", "1620v8", "1622v8", "1623v8", "1625v8", "1626v8", "1630v8",
            "1632v8", "1633v8", "1634v8", "1635v8", "1636v8", "1637v8", "1639v8",
            "1641v8", "1645v8", "1646v8", "1651v8", "1652v8", "1653v8", "1654v8",
            "1657v8", "1658v8", "1659v8", "1663v8"],
        "SP(n = 30)": ["ES346", "ES350", "ES352", "ES353", "ES357", "ES360", "ES363",
            "ES365", "ES366", "ES373", "ES374", "ES389", "ES394", "ES406",
            "ES411", "ES414", "ES415", "ES418", "ES420", "ES424", "ES425",
            "ES426", "ES432", "ES437", "ES457", "ES459", "ES494", "ES499", "ES501"],
        "SP(n = 20)_HX": ["ES344", "ES356", "ES371", "ES372", "ES375", "ES379", "ES382",
            "ES393", "ES396", "ES433", "ES438", "ES440", "ES450", "ES451",
            "ES452", "ES455", "ES487", "ES488", "ES489", "ES556", "ES587"]
    }

    reordered_columns = []
    for group in sample_groups.values():
        reordered_columns.extend([col for col in group if col in scaled_heatmap_data.columns])
    scaled_heatmap_data = scaled_heatmap_data[reordered_columns]

    # Calculate dynamic figure height based on number of peptides
    height = min(max(4, len(heatmap_data) * 0.3), 10)
    plt.figure(figsize=(10, height))

    # Create the heatmap
    heatmap = sns.heatmap(
        scaled_heatmap_data,
        cmap="Blues",
        cbar_kws={'label': 'Binary positivity (0 or 1)'},
        xticklabels=False,
        yticklabels=True,
    )

    current_position = 0
    group_boundaries = []
    for samples in sample_groups.values():
        size = len([s for s in samples if s in scaled_heatmap_data.columns])
        if size > 0:
            current_position += size
            group_boundaries.append(current_position)

    for boundary in group_boundaries[:-1]:
        plt.axvline(x=boundary, color='red', linestyle='--', linewidth=1)

    # Calculate automatic offset based on the number of peptides in the heatmap
    # The more peptides, the more negative the offset needs to be
    if label_y_offset is None:
        num_peptides = len(heatmap_data)
        # Base offset with scaling factor. More peptides = more negative offset
        label_y_offset = -0.6 * (num_peptides / 10) - 0.5
        # Ensure a reasonable minimum offset
        label_y_offset = min(label_y_offset, -0.6)
        print(f"📏 Calculated automatic y-offset: {label_y_offset:.2f}")

    plt.subplots_adjust(top=0.85)
    midpoints = [(group_boundaries[i - 1] + b) // 2 if i > 0 else b // 2 for i, b in enumerate(group_boundaries)]
    
    # Add group labels at the automatically calculated offset
    for name, midpoint in zip(sample_groups.keys(), midpoints):
        plt.text(midpoint, label_y_offset, name, ha='center', va='center', fontsize=10, color='black', rotation=90)

    # Remove title and just use y-axis label for regions
    plt.ylabel(f"{protein_name} Regions (Start_End)")
    
    # Create a combined label with "Samples" and the protein name
    plt.xlabel(f"Samples\n{protein_name}", fontweight='bold')
    
    # Adjust the figure layout to accommodate the expanded xlabel
    plt.tight_layout()
    
    # Create output directory if it doesn't exist
    output_dir = "protein_heatmaps"
    os.makedirs(output_dir, exist_ok=True)
    
    # Use sanitized filename for saving
    safe_filename = sanitize_filename(protein_name)
    output_path = os.path.join(output_dir, f"{safe_filename}_heatmap.png")
    
    plt.savefig(output_path, dpi=300, bbox_inches='tight')
    print(f"✅ Saved heatmap to: {output_path}")
    plt.show()

# Load data
file_path = r"C:\Users\gantak\OneDrive - UMass Chan Medical School\Documents\virscan Krishna\new data analysis\jan 2025\filtered_annotated_consolidated_protein_scores_2_5_both_av.csv"
data = pd.read_csv(file_path)

# Protein lengths
proteins = {
    "BALF2": 1128,
    "BGLF3.5": 153,
    "BZLF1": 245,
    "Capsid protein P40": 605,
    "Capsid protein VP26": 176,
    "DNA replication helicase (EC 3.6.4.-)": 809,
    "EBNA-1": 641,
    "EBNA-6": 1069,
    "EBNA-LP": 506,
    "EBNA2": 487,
    "EBNA3": 944,
    "EBNA3B": 915,
    "EBNA_4": 938,
    "Early antigen protein D": 404,
    "Glycoprotein 42 (gp42)": 223,
    "Gp350": 907,
    "Immediate-early protein Rta": 605,
    "LMP1": 497,
    "Large tegument protein (LTP)": 3154,
    "Major tegument protein (MTP) (Protein p140)": 1318,
    "Major_capsid_protein ": 1381,
    "Protein BDLF2": 420,
    "Protein BMRF2": 357,
    "Protein BOLF1": 1240,
    "Protein LF2": 429,
    "Secreted protein BARF1 (33 kDa early protein) (p33)": 221,
    "Serine/threonine-protein kinase BGLF4 (EC 2.7.11.1)": 429,
    "Shutoff alkaline exonuclease (EC 3.1.11.-)": 470,
    "Tegument protein BKRF4": 217,
    "Tegument protein BRRF2": 537,
    "Tegument_protein_BLRF2": 162,
    "Tegument_protein_BSRF1": 218,
    "Thymidine kinase (EC 2.7.1.21)": 607,
    "Transcriptional activator BRRF1": 310,
    "Triplex capsid protein VP19C homolog": 364,
    "Triplex capsid protein VP23 homolog": 301,
    "Uracil-DNA glycosylase (UDG) (EC 3.2.2.27)": 255,
    "Viral interleukin-10 homolog (vIL-10) (20 kDa protein) (Protein BCRF1)": 170,
    "Virion egress protein BFLF2 (Primary envelopment factor BFLF2)": 318,
    "Virion egress protein UL34 homolog (Primary envelopment factor UL34 homolog)": 336,
    "Virion egress protein UL7 homolog": 278,
    "Virion-packaging protein UL17 homolog": 507,
    "Virion-packaging protein UL25 homolog": 570,
    "gB": 857,
    "gH": 706,
    "gM": 405
}

label_offsets = {
    'BALF2': -5.77,
    'BGLF3.5': -3.6,
    'BZLF1': -4.37,
    'Capsid protein P40': -5.0,
    'Capsid protein VP26': -3.0,
    'DNA replication helicase (EC 3.6.4.-)': -4.85,
    'EBNA-1': -4.0,
    'EBNA-6': -6.0,
    'EBNA-LP': -3.0,
    'EBNA2': -4.72,
    'EBNA3': -5.66,
    'EBNA3B': -5.49,
    'EBNA_4': -5.63,
    'Early antigen protein D': -4.0,
    'Glycoprotein 42 (gp42)': -3.4,
    'Gp350': -5.44,
    'Immediate-early protein Rta': -4.2,
    'LMP1': -3.98,
    'Large tegument protein (LTP)': -18.92,
    'Major tegument protein (MTP) (Protein p140)': -7.91,
    'Major_capsid_protein ': -8.29,
    'Protein BDLF2': -4.52,
    'Protein BMRF2': -4.14,
    'Protein BOLF1': -7.44,
    'Protein LF2': -4.57,
    'Secreted protein BARF1 (33 kDa early protein) (p33)': -2.33,
    'Serine/threonine-protein kinase BGLF4 (EC 2.7.11.1)': -3.57,
    'Shutoff alkaline exonuclease (EC 3.1.11.-)': -2.82,
    'Tegument protein BKRF4': -2.3,
    'Tegument protein BRRF2': -4.22,
    'Tegument_protein_BLRF2': -1.97,
    'Tegument_protein_BSRF1': -2.31,
    'Thymidine kinase (EC 2.7.1.21)': -4.0,
    'Transcriptional activator BRRF1': -2.86,
    'Triplex capsid protein VP19C homolog': -3.18,
    'Triplex capsid protein VP23 homolog': -2.81,
    'Uracil-DNA glycosylase (UDG) (EC 3.2.2.27)': -2.53,
    'Viral interleukin-10 homolog (vIL-10) (20 kDa protein) (Protein BCRF1)': -1.62,
    'Virion egress protein BFLF2 (Primary envelopment factor BFLF2)': -2.91,
    'Virion egress protein UL34 homolog (Primary envelopment factor UL34 homolog)': -3.02,
    'Virion egress protein UL7 homolog': -2.67,
    'Virion-packaging protein UL17 homolog': -3.44,
    'Virion-packaging protein UL25 homolog': -4.42,
    'gB': -5.14,
    'gH': -4.64,
    'gM': -3.43
}

# Generate heatmaps
for protein_name, protein_length in proteins.items():
    print(f"\nGenerating heatmap for {protein_name} (length: {protein_length})")
    # Pass None for automatic offset calculation instead of using the predefined offsets
    generate_protein_heatmap(data, protein_name, protein_length, label_y_offset=None)